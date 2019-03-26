# Exit if any command fails (since these only make sense in sequence, we'd get unpredictable behavior otherwise)
set -e

# Stash any current changes.
git stash

# Validate tag
if git tag | grep -q $1; then
  echo "You passed in a valid tag."
  # Checkout that tag (ensures that the right docs directory corresponds to the right tag)
  git fetch
  git checkout tags/$1
  echo "Checked out tag $1."
  target=$1
elif [ "$1" = "master" ]; then
  echo "You selected the master branch."
  git checkout master
  git pull
  target="dev"
else
  echo "Sorry, you didn't pass in a valid tag."
  exit 1
fi

# Weave examples
(cd docs/examples && julia generate.jl) # weaves notebooks to docs/examples. requires weave to be installed
echo "Weaved example notebooks."

# Run pdflatex steps
(cd docs/tex && pdflatex discretized-differential-operator-derivation.tex && pdflatex discretized-differential-operator-derivation.tex) # compiled pdf to docs/tex
echo "Compiled PDF."

# Clone gh-pages branch
git clone -b gh-pages --single-branch https://github.com/QuantEcon/SimpleDifferentialOperators.jl
echo "Cloned docs branch."

# Pass files to docs repo
mkdir -p SimpleDifferentialOperators.jl/$target/generated/
mv docs/tex/discretized-differential-operator-derivation.pdf SimpleDifferentialOperators.jl/$target/generated/ # move PDF
echo "Copied PDF."
mv docs/examples/*.ipynb SimpleDifferentialOperators.jl/$target/generated/ # move notebooks
echo "Copied example notebooks."
mv docs/examples/*.html SimpleDifferentialOperators.jl/$target/generated/ # move notebooks
echo "Copied example HTML."

# Git operations
(cd SimpleDifferentialOperators.jl; git add -A; git commit -m "Add generated objects to $target docs"; git push)
echo "Carried out git operations."

# Clean
rm -rf SimpleDifferentialOperators.jl
git checkout **/*.pdf
echo "Cleaned."
