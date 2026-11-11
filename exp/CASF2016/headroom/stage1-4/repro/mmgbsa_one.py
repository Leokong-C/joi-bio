
import sys, os, json, time, pathlib
import numpy as np
code = sys.argv[1]
rec = {"code": code}
try:
    from openff.nagl_models import _dynamic_fetch as DF
    pathlib.Path(DF.CACHE_DIR).mkdir(parents=True, exist_ok=True)
    import openmm, openmm.app as app, openmm.unit as unit
    from pdbfixer import PDBFixer
    from openff.toolkit import ForceField, Topology, Molecule
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
    from openff.units import unit as off_unit

    OBC = {1:(0.12,0.85), 6:(0.17,0.72), 7:(0.155,0.79), 8:(0.15,0.85), 9:(0.15,0.88),
           15:(0.185,0.86), 16:(0.18,0.96), 17:(0.17,0.80), 35:(0.20,0.80), 53:(0.20,0.80)}
    LIGPATH = {"3syr": ("sdf","casf/lig_fix/3syr.sdf"), "1a30": ("pdb","casf/lig_fix/1a30_pep.pdb")}

    def add_gb(system, Z):
        nb = [system.getForce(i) for i in range(system.getNumForces())
              if isinstance(system.getForce(i), openmm.NonbondedForce)][0]
        gb = openmm.GBSAOBCForce()
        gb.setSolventDielectric(78.5); gb.setSoluteDielectric(1.0)
        gb.setNonbondedMethod(openmm.GBSAOBCForce.NoCutoff)
        gb.setSurfaceAreaEnergy(2.25936*unit.kilojoule_per_mole/unit.nanometer**2)
        for i, z in enumerate(Z):
            q = nb.getParticleParameters(i)[0]
            r, s = OBC.get(z, (0.15, 0.80))
            gb.addParticle(q, r*unit.nanometer, s)
        system.addForce(gb)
        return system

    def energy(system, Z, positions, minimize=0):
        system = add_gb(system, Z)
        integ = openmm.VerletIntegrator(1*unit.femtosecond)
        ctx = openmm.Context(system, integ, openmm.Platform.getPlatformByName("CPU"))
        ctx.setPositions(positions)
        if minimize:
            openmm.LocalEnergyMinimizer.minimize(ctx, maxIterations=minimize)
        st = ctx.getState(getEnergy=True, getPositions=True)
        e = st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
        pos = st.getPositions(asNumpy=True)
        del ctx, integ
        return e, pos

    t0 = time.time()
    fixed = f"mmgbsa/fixed/{code}.pdb"
    if not os.path.exists(fixed):
        fx = PDBFixer(filename=f"casf/prot/{code}.pdb")
        fx.findMissingResidues(); fx.missingResidues = {}
        fx.findNonstandardResidues(); fx.replaceNonstandardResidues()
        fx.removeHeterogens(keepWater=False)
        fx.findMissingAtoms(); fx.addMissingAtoms(); fx.addMissingHydrogens(7.4)
        app.PDBFile.writeFile(fx.topology, fx.positions, open(fixed, "w"), keepIds=True)

    ptop = Topology.from_pdb(fixed)
    allmols = list(ptop.molecules)
    # PDBFixer's removeHeterogens keeps free standard residues that sit in ATOM records as
    # their own chain (e.g. the zwitterionic proline in 3bv9) — crystallisation additives, not
    # receptor. ff14sb has no template for an isolated zwitterion, so openff falls through to
    # the small-molecule handler and raises NonIntegralMoleculeChargeError. Drop them.
    pmols = [m for m in allmols if m.n_atoms >= 60]
    dropped = [{"n_atoms": int(m.n_atoms), "smiles": m.to_smiles(explicit_hydrogens=False)}
               for m in allmols if m.n_atoms < 60]
    if not pmols:
        raise RuntimeError("no receptor chain above 60 atoms")
    nP = sum(m.n_atoms for m in pmols)

    fmt, path = LIGPATH.get(code, ("sdf", f"casf/lig_sdf/{code}.sdf"))
    # Several CASF ligands have unspecified stereocentres in the deposited SDF (1lpg, 1sqa).
    # MM parameters do not depend on the stereo assignment, so proceed rather than fail.
    lg = Molecule.from_file(path, file_format="PDB" if fmt == "pdb" else "SDF",
                            allow_undefined_stereo=True)
    if isinstance(lg, list):
        lg = lg[0]
    NAGLToolkitWrapper().assign_partial_charges(lg, partial_charge_method="openff-gnn-am1bcc-1.0.0.pt")
    # NAGL charges do not sum exactly to the formal charge; spread the residual evenly,
    # otherwise create_interchange raises NonIntegralMoleculeChargeError (seen on zwitterions).
    q = lg.partial_charges.m_as(off_unit.elementary_charge)
    resid = float(lg.total_charge.m) - float(q.sum())
    lg.partial_charges = (q + resid/len(q)) * off_unit.elementary_charge

    FF = ForceField("ff14sb_off_impropers_0.0.4.offxml", "openff-2.2.0.offxml")
    comb = Topology.from_molecules(pmols + [lg])
    assert comb.n_atoms == nP + lg.n_atoms
    iC = FF.create_interchange(comb, charge_from_molecules=[lg])
    sC = iC.to_openmm(combine_nonbonded_forces=True)
    sP = FF.create_interchange(Topology.from_molecules(pmols)).to_openmm(combine_nonbonded_forces=True)
    sL = FF.create_interchange(Topology.from_molecules([lg]),
                               charge_from_molecules=[lg]).to_openmm(combine_nonbonded_forces=True)
    Z = [a.atomic_number for a in comb.atoms]
    pos = np.array(iC.positions.m_as("nanometer"))
    eC, pmin = energy(sC, Z, pos*unit.nanometer, minimize=200)
    pm = pmin.value_in_unit(unit.nanometer)
    eP, _ = energy(sP, Z[:nP], pm[:nP]*unit.nanometer, 0)
    eL, _ = energy(sL, Z[nP:], pm[nP:]*unit.nanometer, 0)
    rec.update(E_complex=round(eC, 3), E_receptor=round(eP, 3), E_ligand=round(eL, 3),
               dG_mmgbsa=round(eC-eP-eL, 3), n_prot=nP, n_lig=int(lg.n_atoms),
               n_chains=len(pmols), n_dropped=len(dropped), dropped=dropped, net_charge_residual=round(resid, 4),
               sec=round(time.time()-t0, 1))
except Exception as exc:
    rec["error"] = f"{type(exc).__name__}: {exc}"[:200]
print("@@RESULT@@" + json.dumps(rec))
