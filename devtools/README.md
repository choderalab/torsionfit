Developer Notes / Tools / License
=================================

Assorted notes for developers.

How to do a release
-------------------
- Update the `docs/whatsnew.rst` document. Use the github view that shows all the
  commits to master since the last release to write it.
- Update the version number in `setup.py`, change `ISRELEASED` to `True`
- Update the version number in `devtools/conda-recipe/meta.yaml`
- Commit to master, and [tag](https://github.com/choderalab/torsionfit/releases) the
  release on github
- To push the source to PyPI, use `python setup.py sdist --formats=gztar,zip upload`.
  **Before pushing to PyPI, run `git clean -x -f -d` to make sure that nothing
  else gets added to the source release by accident**
- Update `meta.yaml` in the omnia-md/conda-recipes repository. Commit
  to master. Travis and Appveyor will build binary conda packages.
- After tagging the release, make a NEW commit that changes `ISRELEASED` back
  to `False` in `setup.py`


It's important that the version which is tagged on github for the release be
the one with the ISRELEASED flag in setup.py set to true.

