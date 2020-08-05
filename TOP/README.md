# CMS-TOP-Rivet

This GitLab repository contains Rivet routines for Top PAG related analyses and configuration files to run them. 
Routines are gather in the `CMS-TOP-Rivet/RivetTop` subpackage. Some may not be validated.

## Organisation of `Rivet/Top` 

The `master` branch contains all the plugins and data.
The `scr/` and `interface/` repositories contain `.cc` and `.h` files respectively for these plugins.
The `data/` repository contains `.yoda` (for data and MC) and `.plot` files.
The `test/` repository contains only small scripts that run on `.yoda` files (for corrections, normalization...).

Common interfaces to be used in analysis modules are available.

  * PartonTop : Build parton level top quarks and its decay products.
  * PseudoTop : Build particle level top quarks (pseudo top) and its decay products.
  * CMSGenParticle (not validated) : Emulate CMS genParticlesForJets in Run-I

Some analyses done in CMS are also available.

  * TOP-12-028 (CMS\_2015\_I1370682.cc) : Differential cross section measurement, under validation
  * TOP-12-042 (CMS\_TOP\_12\_042.cc) : Differential cross section using event variables, under validation

Analysis plugins that rely on parton level information are developed for internal use and correction function extraction.
  * TOP-12-028 (CMS\_2015\_I1370682\_internal.cc) : Differential cross section (synchronized with MadGraph in the paper)
  * TOP-12-041 (CMS\_TOP\_12\_041\_internal.cc) : Jet multiplicity measurement (under validation)

Other modules are also available for the MC studies

  * Les Houches working group (CMS\_LesHouches2015.cc) : Analysis for the Les Houches ttbar working group using the dilepton channel (see separate [documentation page](doc/LesHouches.md))
  * MC\_TTBAR\_HADRON.cc : Lepton+jets channel study
