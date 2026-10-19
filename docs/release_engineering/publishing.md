# Publishing and post-release verification

`publish.yml` is the only supported PyPI promotion route. It is manual and runs
from the default branch. It does not build new distributions. It downloads the
CI distribution artifact named by the committed version-specific approval
record, verifies the CI and Scientific Tier 2 run identity and success, checks
both artifact hashes, and carries the exact files through a retained promotion
bundle. The protected `pypi` environment is required before the tag is created
or the PyPI OIDC publisher runs.

## Before dispatch

1. Merge the candidate commit into the default branch without rewriting the
   candidate commit. Record the exact candidate SHA, version, successful CI and
   Tier 2 run IDs, expected job inventory, distribution filenames and hashes,
   scientific scope and exclusions, and owner approval in
   `docs/release_engineering/releases/<version>.json`.
2. Confirm the version-specific PyPI endpoint returns 404 and that the version
   tag does not exist. Never reuse or move a release tag.
3. Dispatch `publish.yml` on the default branch with the version, candidate
   SHA and exact successful run IDs from the approval record.
4. Review and approve the pending `pypi` environment deployment. The workflow
   rechecks all approval fields and artifact hashes inside the protected job,
   then creates `v<version>` at the candidate SHA and uploads those verified
   bytes through PyPI Trusted Publishing.

The dispatch job fails closed for a different ref, candidate, workflow run,
job set, version, artifact set, artifact hash, existing PyPI version, or
existing tag. Tag pushes do not invoke a second publisher.

## Post-publication

The final workflow job polls the version-specific PyPI JSON endpoint, verifies
that the wheel and sdist names and SHA-256 hashes match the approval record,
installs that exact version in a fresh virtual environment, runs `pip check`,
checks the installed import origin and version, and runs `pyar-cli --help` and
`pyar --help`. A passing upload job alone is not release completion; G7 requires
the post-publication job to pass.

If upload succeeds but G7 fails, do not move the tag or upload replacement
files. Preserve the failed run, artifact bundle and registry response, notify
the project owner, and decide whether to document the incident, yanking, or a
new patch version. PyPI artifacts are immutable after upload.

The environment's required reviewer and deployment branch restriction are
repository settings and must be verified independently of this workflow file.
The repository's release record must retain the workflow URL, tag commit, PyPI
file hashes and clean-install result.
