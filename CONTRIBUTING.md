# Basic workflow

* Fork the repository
* Add a preliminary version of your plugin
* Create a work-in-progress (WIP) merge request
* Receive comments and finalize plugin
* Add validation plots to the merge request using the [cmsRivetPlots2Web](https://gitlab.cern.ch/cms-gen/Rivet/blob/master/scripts/cmsRivetPlots2Web) script
* If code and validation is ok, CMS Rivet contacts will send the plugin to the Rivet authors

#### Resolving merge conflicts

If your merge request shows a message `This merge request contains merge conflicts`, please follow these steps:

1. Add the cms-gen main repository to your .git/config:
```
[remote "cms-gen"]
        url = https://:@gitlab.cern.ch:8443/cms-gen/Rivet.git
        fetch = +refs/heads/*:refs/remotes/cms-gen/*
```
2. `git fetch cms-gen`
3. `git merge cms-gen/master`
4. After resolving the conflicts, commit and push again.

This can also be used to sync your fork with the cms-gen repository anytime.

# Naming convention

* Papers: `CMS_{YYYY}_I{inspire-id}.cc`
* PAS: `CMS_{YYYY}_PAS_TOP_{YY}_{NNN}.cc`

If you do not have an inspire-id for your paper yet, please use the PAS naming convention.

If your PAS has an inspire-id, don't use it. Use the PAS naming convention.

# Plot style

Add `Title=CMS, xx TeV, process [<extra information>]` to your plots to make the CMS origin easily identifiable if theorists use the plugin.

# Coding style

Please follow the Rivet coding style documented here:
* https://rivet.hepforge.org/trac/wiki/WritingAnAnalysis
* https://rivet.hepforge.org/trac/wiki/CodingStyleGuide
