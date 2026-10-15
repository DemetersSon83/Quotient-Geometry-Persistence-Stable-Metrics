# GitHub setup

The supplied archive contains an initialized Git repository on branch `main`, with one packaging commit and no remote. The original calculation scripts and saved study results are included; Rocq files are excluded.

Create an empty repository under the intended GitHub account or organization. Choose its visibility explicitly. Leave GitHub's README, license and gitignore initialization options unchecked so that the empty remote can receive this repository's existing commit.

Extract the archive, open a terminal inside `quotient-geometry-reproducibility`, and run:

```sh
git status
python3 reproduce.py
git remote add origin https://github.com/OWNER/REPOSITORY.git
git push -u origin main
```

Replace `OWNER/REPOSITORY` with the real destination. Authentication uses your normal GitHub credentials; no credential is included in this package. On Windows, `python` may be the command for your Python 3.12 interpreter.

The workflow in `.github/workflows/reproduce.yml` runs full reproduction on pushes, pull requests, and manual dispatch. It has read-only repository permissions and no external Python dependencies. Remote workflow execution has not yet been tested; the committed `validation/` files document the local run.

The existing study archive did not specify a repository-wide license. `NOTICE.txt` preserves that status; the rights holder can select a license when deciding release terms. Update repository URLs or archival identifiers in `CITATION.cff` once they exist.

If using a source-only export that omits Git history, first initialize and commit the files with `git init -b main`, `git add .`, and `git commit -m "Package manuscript reproduction code"` before adding the remote.
