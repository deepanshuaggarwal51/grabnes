# Jeil Jung Group Codes

Let's organize all our codes and notebooks into a common git folder. 

Let me know if you want to add other folders. Let's try to keep it organized.

If you want to do code development branching off the main branch (highly recommended rather than working in your own untracked version), ask me to make you a branch in your name. Regularly commiting to github will also increase your market value if you ever need to apply for an IT job in industry...

Note, only upload codes that are ok to share with everyone. If it's a code from a collaborator, make sure he is ok to share it here.

Folders (for now):
- Matlab codes from Jeil.
- Lanczos code from Nicolas (based on open source version by Rafael Martinez Gordillo)
- PyBinding related notebooks. This includes the projects investigated by the students during the summer as well as the final version of the students having worked on it during the visual physics class. Refer to the official github for the PyBinding source code. Note that the version installed by Pip is not the last PyBinding version. To run the twistedBilayer graphene notebook for instance, you need the github version.
- Small scripts that can be useful for several projects.

Other possible folders
- LAMMPS input scripts and readme files by Bheema
- Useful .bashrc, .vimprofile, etc type of files
- Makefiles for different environments
- Machine learning folder

How to clone and work with this folder:
- click on the clone folder button of this repository and copy the address
- type "git clone addressOfRepository" in your terminal where you want to create this repository
- if you have your own branch (and you should), make sure to type "git checkout nameOfBranch" the first type you start using it (otherwise, you will make changes in the master brach, which might create version conflicts with other users).
- type "git pull" at the beginning of EACH time you want to work in the rep
- do what you want in your local directory
- type the following commands at the end of each work session, in order: "git add nameOfFileYouHaveCreatedOrModified", git commit -m "a descriptive message of the changes you have made", "git push"

Other useful git commands:
- "git status" to get useful information during from your git work session
- "git diff nameOfFiles" to see how files have changed


Note from github rules:

We recommend repositories be kept under 1GB each. This limit is easy to stay within if large files are kept out of the repository. If your repository exceeds 1GB, you might receive a polite email from GitHub Support requesting that you reduce the size of the repository to bring it back down.

In addition, we place a strict limit of files exceeding 100 MB in size. For more information, see "Working with large files."

Based on this, I would suggest only adding code and notebooks. Keep the number of figure output to a minimum. Large data files should also be avoided. One can for instance use os.chdir("/directory/to/readOrWrite/data/fromTo/") in python to link your notebooks to outside folders in either our shared Dropbox or your personal folders.
