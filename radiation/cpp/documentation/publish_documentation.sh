#!/bin/bash -e

# Quit on errors.
set -o errexit -o nounset

# Settings.
DOCUMENTATION_PATH=radiation/cpp/documentation
CHANGESET=$(git rev-parse --verify HEAD)

# Make documentation.
doxygen Doxyfile

# Make sure branches are up to date.
git remote set-branches --add origin gh-pages
git fetch origin

# Commit documentation in master.
cd ../../..
git add ${DOCUMENTATION_PATH}
git commit -m "adding documentation"

# Check out gh-pages branch and merge documentation from master commit.
git checkout gh-pages
rm -r ${DOCUMENTATION_PATH}
git checkout master ${DOCUMENTATION_PATH}

# Add the merged changes and push.
git commit -a -m "Automated documentation build for changeset ${CHANGESET}."
git push -u origin gh-pages

# Checkout master again and blow away generated docs.
git checkout master
git rm ${DOCUMENTATION_PATH}/html/ ${DOCUMENTATION_PATH}/latex/

echo "-- Successfully updated documentation!"
