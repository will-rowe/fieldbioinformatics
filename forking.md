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

Pull the desired branch from the upstream repository. This method will retain the commit history without modification.
```bash
git pull -m "merging" https://github.com/artic-network/fieldbioinformatics.git master
```

Push the merge to your GitHub repository
```bash
git push origin master
```
