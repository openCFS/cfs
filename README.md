CFS++
==============================

This is the git repository of CFS++.
It is currently kept in sync with our [SVN repo](https://cfs.mdmt.tuwien.ac.at/svn/CFS++).
This is done by [subgit](https://subgit.com/).
Integration ws done accoding to the [gitlab-guide](https://docs.gitlab.com/ce/user/project/import/svn.html).
Our subgit instance is registered with the following data, which was issued for free (academic project) by TMate Software
```
Registered for:       CFS++ (http://cfs-doc.mdmt.tuwien.ac.at)
Purchase ID:          OS-515510923183017
Expiration date:      November 17, 2018
```

It should be used for testing, i.e. for writing developer docu in the branch `devel_docu`.

[Developer Documentation](/share/doc/developer/README.md)
---------------

Working locally
---------------

1. clone the repository `git clone ...`
2. checkout the branch `git checkout devel_docu`
3. add, commit & push

To avoid unnecessary merge commits, tell git to always rebase when pulling, to do this on a project level add this to your `.git/config` file
```
[branch "devel_docu"]
  rebase = true
```
Of course there are [different possibilities to avoid merge commits](http://kernowsoul.com/blog/2012/06/20/4-ways-to-avoid-merge-commits-in-git/).

Working on gitlab
-----------------

You can edit or create files directly on gitlab.
This way one can use the `preview` function of the editor.
Be sure to
* write commit messages
* commit to the correct branch 
