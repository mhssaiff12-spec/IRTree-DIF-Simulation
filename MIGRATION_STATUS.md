# IRTree-DIF: migration status

The connected ChatGPT session established write access to this repository on 2026-10-16.

## Current status
- Repository exists. **The 34 MB portable simulation package has NOT been uploaded to GitHub.**
- A manual-only verification workflow is installed: `.github/workflows/verify-irtree.yml`.
- No full Monte Carlo simulations have been launched in GitHub Actions.
- Existing Claude Cloud analysis must not be stopped until data are backed up and verified.
- The corrected source package is `IRTree_DIF_portable_core8_GitHub_safe.zip`, available as a downloadable artifact in the ChatGPT conversation.

## One-time upload requirement
1. Download the safe ZIP on the user's computer.
2. Extract the ZIP and upload the *contents of the portable folder*, keeping folders like `sim/`, `R_env/`, `scripts/`, and `.github/` at repository root.
3. Do not upload personal information or secrets; this repository is public.
4. Avoid overwriting `.github/workflows/verify-irtree.yml`. Do not enable the long-running workflow until the verification passes.
5. Before running, obtain the latest Claude results backup, not just the nine-replication snapshot.

## Verification
Go to Actions > IRTree verification (manual only) > Run workflow. Check for `PASS` before enabling any simulation.

## Precautions
- Pin R 4.3.3 and package versions from the source bundle.
- Linux versus Claude numerical reproducibility needs actual testing.
- GitHub artifact storage and branch persistence are independent; do not rely on artifacts as permanent backups.
- The main workload uses substantial CPU resources and must be explicitly enabled after verification and backup.
