1. Create new empty repo in GitHub (say, `artic-ebov`)

2. Clone this repo locally:
```bash
git clone git@github.com:artic-network/fieldbioinformatics
```

3. Push it to new repo from working copy of original:
```bash
cd fieldbioinformatics
git remote set-url origin git@github.com:artic-network/artic-ebov
git push origin master
```

4. Clone a local copy of the new repo:
```bash
git clone git@github.com:artic-network/artic-ebov
```

### Merging with the original repo

Add fieldbioinformatics as the upstream repo:
```bash
cd artic-ebov
git remote add upstream https://github.com/artic-network/fieldbioinformatics
```

Fetch the branches and their respective commits from the upstream repository. Commits to master will be stored in a local branch, upstream/master.
```bash
git checkout master
git merge upstream/master
```

Merge the changes from upstream/master into your local master branch. This brings your fork's master branch into sync with the upstream repository, without losing your local changes.
```bash
git merge upstream/master
```
