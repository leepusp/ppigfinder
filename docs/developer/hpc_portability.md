# Portable AlphaFold 3 / HPC integration contract

**Status: implementation contract and audit, not a claim of completed automation.** Reviewed against source files on 2026-11-18.

## Goal

A ppigFinder installation must be useful without DaVinci, Slurm, SSH or local GPUs. Cluster functionality is an optional adapter. The core application owns selection, validated job specifications, manifests, outputs and result parsing; the HPC backend owns transport and scheduling. AlphaFold 3 code, model weights and databases are **not bundled**.

## Current code map

| Responsibility | Current code | Verified code capability |
| --- | --- | --- |
| CLI wrapper command | `ppigfinder/alphafold/af3_cli.py` | Builds command argv and probes `af3 --help` |
| Direct AF3 service | `ppigfinder/services/alphafold_service.py` | Previews/exports/runs a local `af3` command |
| AF3 input writer | `ppigfinder/structure_prediction/af3_input_writer.py` | Writes AF3 JSON, initially protein-only |
| Batch orchestration | `ppigfinder/services/structure_prediction_service.py` | Plans jobs, creates artifacts/scripts |
| Slurm scripts | `ppigfinder/structure_prediction/slurm_batch_writer.py` | Arrays and task tables; runs wrapper with `--executor local` |
| Cluster resource hints | `ppigfinder/structure_prediction/cluster_profiles.py` | `generic_slurm` and optional `davinci` mappings |
| Result parsing | `ppigfinder/services/alphafold_results_service.py` | Sequential/parallel folder parsing |
| Modular SSH/monitoring | `ppigfinder/hpc/ssh_client.py`, `job_monitor.py`, `ppigfinder/services/hpc_service.py` | Empty in reviewed revision; do not advertise complete automation |
| Stable desktop | `ppigfinder/app.py` | Loads active `legacy_v29_14.py` GUI |

## Existing execution patterns to keep distinct

1. **Wrapper-managed Slurm:** `AlphaFoldService.build_slurm_options()` produces `af3` flags such as `--slurm-partition`, `--slurm-gres`, etc. This mode allows the wrapper to submit/manage resources when supported.
2. **Outer Slurm allocation:** `slurm_batch_writer.py` creates an `sbatch` script. Inside that allocated job, it calls `af3 --json-path ... --stage all --executor local`. **Never insert wrapper Slurm submission flags here** unless the wrapper explicitly requires them.
3. **Manual server export:** AlphaFold Server JSON is a separate dialect/workflow and must not be blindly passed to a local runner.

## Environment contract (portable runner)

Each institution must document these site-specific settings:

| Field | Purpose | Default policy |
| --- | --- | --- |
| `name` | Profile identifier | `generic_slurm` |
| `ssh_host`, `ssh_port`, `ssh_user` | Connection | Must be supplied for remote transport |
| `jump_host` | Optional bastion | Disabled by default |
| `remote_root` | Cluster-visible workspace | Required for remote submission |
| `scheduler` | Scheduler adapter | `slurm`, otherwise manual |
| `partition`, `account`, `qos` | Scheduling policy | Site-defined; never hard-code DaVinci |
| `cpus`, `memory`, `walltime` | Resource plan | Estimate, then review per site |
| `gpu_request` | Actual scheduler directive | Site-defined; `gres/shard` is not universal |
| `environment_setup` | Module/conda/container | User/admin approved commands |
| `af3_command` | Executable or installed wrapper | `af3` only if on PATH |
| `model_dir`, `db_dir`, `image` | Backend-specific inputs | Site-defined; avoid private paths in shared docs |
| `stage`/`executor` | Invocation semantics | Match the detected wrapper |

Do not store secrets, private keys, passwords, database credentials or personal home paths in committed cluster profiles. Prefer SSH agents/config and local nontracked configuration files.

## Technical defects / caveats found during review

- `slurm_batch_writer.py` currently defaults to `module load alphafold3/1.0`, which is **not portable** and contradicts `ppigfinder/hpc/profiles/davinci.py`, which assumes `af3` is already on PATH. The provided `scripts/hpc/af3_slurm_generic.sh` is a **reference template**, not a connected fix in the existing GUI.
- Generated task tables contain **absolute paths**. Copying them from a laptop to an HPC server without rewriting/re-rendering paths will cause failures. Plan artifacts directly under the remote root or provide an explicit path-rebasing step.
- `AF3_SCRIPT` is a configurable command in Slurm scripts, but module/environment initialization is currently hard-coded unless the caller passes an explicit `module_commands` override.
- The resource planner maps generic GPU size classes to DaVinci-specific `shard:10/32/50`; these values must be checked against the current Slurm configuration before use.
- The package's `hpc_service.py`, `ssh_client.py` and `job_monitor.py` do not implement the promised remote lifecycle. Legacy GUI may use other routines; those should be audited before migration.
- Current AF3 writer supports proteins, not complete AF3 DNA/RNA/ligand functionality.

## Suggested interfaces for the next implementation PR

```python
class HPCTransport:
    def test_connection(self): ...
    def upload_directory(self, local_path, remote_path): ...
    def download_directory(self, remote_path, local_path): ...

class SchedulerAdapter:
    def submit(self, remote_script_path): ...
    def status(self, job_id): ...
    def cancel(self, job_id): ...
    def logs(self, job_id): ...

class StructureBackendAdapter:
    def validate_environment(self): ...
    def build_inputs(self, spec): ...
    def render_command(self, spec, execution_context): ...
    def parse_results(self, path): ...
```

Implement separate `LocalTransport`, `SSHTransport`, `SlurmScheduler`, and `ManualScheduler` adapters, with cluster configuration independent of orchestration. Use a stable job manifest recording input SHA-256 checksums, backend version, resource profile, submission timestamp, scheduler job ID, log locations and terminal state.

## Required tests before stating that remote execution is supported

1. Verify install/launch on a clean workstation without `af3`, SSH and Slurm.
2. Validate generated JSON for the specific local AF3 wrapper.
3. Generate one small batch for `generic_slurm` and assert **no site-specific module/partition/GRES** leaks into default output.
4. Use a fake local `af3` executable to test argv and failure propagation without expensive GPU runs.
5. Run `shellcheck` on the generated scripts and check task-table/path quoting.
6. Mock SSH transport and Slurm commands for success, authentication failure, timeout, failed job and resume.
7. Run a controlled real-cluster smoke test with a known job and recorded outputs.
8. Verify result-import round trip and no false biological 'confirmed interaction' status.

Only after those tests should the user manual describe GUI-driven remote submission as automatic.
