## Test environments
* Mac OS Tahoe 26.1, R 4.5.2
* Windows 10 Pro, R 4.5.2
* win-builder (devel and release)
* R-hub: ubuntu-latest, macos-latest, windows-latest on GitHub; 
  macos-15 on GitHub, ASAN + UBSAN on macOS;
  Ubuntu 24.04.3 LTS (-next and -release)
  

## R CMD check results

** Mac OS Tahoe 26.1, R 4.5.2
0 errors | 0 warnings | 0 notes | 1 info
* checking CRAN incoming feasibility ... [2s/12s] INFO
Maintainer: ‘Silvia Castiglione <silvia.castiglione@unina.it>’

Suggests or Enhances not in mainstream repositories:
  rnaturalearthhires
Availability using Additional_repositories specification:
  rnaturalearthhires   yes   http://packages.ropensci.org

-This note is about the package 'rnaturalearthhires' which is available from an additional repository as indicated in the DESCRIPTION of the package.


** Windows 11 Pro, R 4.5.2
0 errors | 0 warnings | 1 note | 1 info
* checking CRAN incoming feasibility ... [10s] INFO
Maintainer: 'Silvia Castiglione <silvia.castiglione@unina.it>'

Suggests or Enhances not in mainstream repositories:
  rnaturalearthhires
Availability using Additional_repositories specification:
  rnaturalearthhires   yes   http://packages.ropensci.org

-This info is about the package 'rnaturalearthhires' which is available from  an additional repository as indicated in the DESCRIPTION of the package.


* checking for future file timestamps ... NOTE
unable to verify current time

-This is an erratic note impossible to resolve locally.


** win-builder (devel and release)
** release and devel
0 errors | 0 warnings | 0 notes | 2 info

* checking CRAN incoming feasibility ... [14s] INFO
Maintainer: 'Silvia Castiglione <silvia.castiglione@unina.it>'

Suggests or Enhances not in mainstream repositories:
  rnaturalearthhires
Availability using Additional_repositories specification:
  rnaturalearthhires   yes   http://packages.ropensci.org


* checking package dependencies ... INFO
Package suggested but not available for checking: 'rnaturalearthhires'

-This infos are about the package 'rnaturalearthhires' which is available from  an additional repository as indicated in the DESCRIPTION of the package.


**  R-hub: ubuntu-latest, macos-latest, windows-latest on GitHub; 
          macos-15 on GitHub, ASAN + UBSAN on macOS; 
          Ubuntu 24.04.3 LTS (-next and -release)
0 errors | 0 warnings | 0 notes 

## Reverse dependencies

The package has no reverse dependency.

