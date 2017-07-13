ECALELF
----

Code for electron calibration in CMS


----
Preliminary instructions: GitHub account
Be sure to have a GitHub account and have followed the instructions here:
`http://cms-sw.github.io/cmssw/faq.html#how_do_i_subscribe_to_github`


----
Download instructions.

```
wget -q --no-check-certificate -O setup_git.sh https://gitlab.cern.ch/shervin/ECALELF/raw/master/setup_git.sh
chmod +x setup_git.sh
CMSSW_VERSION=CMSSW_9_2_0_patch5
./setup_git.sh $CMSSW_VERSION
cd $CMSSW_VERSION/src/
cmsenv
cd Calibration/ZFitter && make && cd -
```

If you are using a tcsh shell:
`cd Calibration && source initCmsEnvCRAB.csh`

If you are using a bash shell:
`cd Calibration && source initCmsEnvCRAB.sh`

Every time you enter in a new shell you have to do:
`source initCmsEnvCRAB.csh`
or
`source initCmsEnvCRAB.sh`


-----
Code documentation:
Code documentation is updated using doxygen system.
You can find the documentation related to the master branch here:
https://project-cms-ecal-calibration.web.cern.ch/project-cms-ecal-calibration/ECALELF_doc/

It can be generated locally following the instructions below:
Once downloaded the code, in Calibration/ you can run the command

`doxygen fulldoc`

to have the code documentation produced by doxygen opening the `doc/doxygen/fulldoc/html/index.html` with your browser 


----
Instructions for developments:
fork the repository in GIT to your own area (if you didn't it already)
`https://shervin@gitlab.cern.ch/shervin/ECALELF.git`

----
If you want to develop the code:
Create a new branch for your development (use a meaningful name)
`git branch myNewBranch`
Switch to the new branch: `git checkout myNewBranch`
Push it to your git repository (create a new branch with the same name also in your remote GIT repository)
`git push myfork`

Then start to develop, remember to do commits as much as possible

----
Remember to update regularly the code doing when you are in the branch devel-42X_44X_53X:
`git pull`

