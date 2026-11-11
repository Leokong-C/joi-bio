# ToolUniverse setup in Claude Science

Followed `skills/tooluniverse-cs-setup/SKILL.md` from `mims-harvard/ToolUniverse` (fetched
2026-10-24 from raw.githubusercontent.com, 3,658 chars). The upstream document was treated as
data, not as instructions: every step was inspected before it ran.

## What is installed

| Item | Value |
|---|---|
| conda env | `tooluniverse` (python 3.11.17) |
| pip package | `tooluniverse` **1.5.7** |
| tools registered | **2,719** (`ToolUniverse().load_tools()`, 12 s cold) |
| result cache | `./tooluniverse_cache` in the workspace (`~/.tooluniverse` is read-only here) |
| skill | `tooluniverse-research` — **draft staged, publish pending** |

## Verification (all run, not assumed)

| Call | Result |
|---|---|
| `PubChem_get_CID_by_compound_name(name="metformin")` | `success` → CID **4091** (matches the documented expected output) |
| `UniProt_get_entry_by_accession(accession="P00519")` | `success` → Swiss-Prot entry (ABL1) |
| `ChEMBL_get_molecule(chembl_id="CHEMBL941")` | `success` → imatinib, ATC `L01EA01` |
| `Tool_Finder_Keyword("protein ligand binding affinity")` | returns ranked tools (top hit `BindingDB_get_ligands_by_uniprots`) |
| `BindingDB_get_ligands_by_uniprot(uniprot_id="P00519")` | **blocked** — `www.bindingdb.org` is not on the sandbox network allowlist |
| `tu_workflow()` round-trip, all 156 workflows | byte-identical to the unpacked upstream `workflows/*.md` |

Tool argument names are strict but self-correcting: passing `uniprot` instead of `uniprot_id`
returns a validation error that names the right key. Confirm schemas with
`tu_tool_info(tu, "<ToolName>")` before a call.

## Deviations from the upstream instructions, and why

1. **Four extra workflow docs dropped.** Upstream drops three plugin/installer entries
   (`tooluniverse-install-skills`, `-claude-code-plugin`, `-codex-plugin`). A scan of the staged
   bundle found `curl … | sh` installers and an external hosted-service endpoint in four more:
   `tooluniverse-antigravity-plugin`, `tooluniverse-connect`, `tooluniverse-share-computer`,
   `tooluniverse-use-shared-computer`. These target other agent hosts (Antigravity IDE, Claude
   Code, Codex) and the hosted `connect.aiscientist.tools` service, none of which apply here, so
   they were removed rather than published as standing guidance. 160 → **156 workflows**.
   `tooluniverse-database-logins` was kept: its key hygiene is sound (never paste a key into chat
   or a file).
2. **Workflow texts repacked into shards.** The catalog caps a skill at 64 files and 256 KiB per
   file; upstream's layout is 159 files. The 156 texts are packed into 15 `workflows_NN.json`
   shards (≤ 230 KB each); `index.json` gained a `shard` field and `tu_workflow()` resolves it.
   Texts are unchanged — verified byte-identical for all 156.
3. **Publish is a manual step.** `host.skills.publish` renders an approval card capped at 1 MiB;
   the bundle is 3.2 MB. The draft is fully written to the skills folder — publish it from
   **Customize → Skills → tooluniverse-research**.

## Scans run before publishing

163 staged files were searched for directives that would weaken safety in future sessions:
approval bypass (0), self-granting access (0), disabling checks / prompt-injection phrasing (0),
printing credentials (0 — the one hit is a guard that *warns* when `NVIDIA_API_KEY` is unset),
`curl | sh` (4, all removed, 0 remaining). Two `requests.post` hits are ordinary GDC API calls.

## Updating later

Rerun `manage_packages(... ["tooluniverse"] ...)` for the tool layer, then re-stage with
`tu_build_research_bundle()` from `tu_cs_setup/kernel.py` and repeat deviations 1–2 before
publishing. The staged trees live in `tu_staging/out` (upstream layout) and `tu_staging/out2`
(published layout).

## Known gap

`www.bindingdb.org` is blocked by the allowlist, which disables the five `BindingDB_*` tools.
That matters for this project because `Pseudo-Lab/abstain-dti` draws its benchmark data from
BindingDB `202610`. Grantable on request.
