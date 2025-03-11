The .gitlab directory
=====================

This directory is read by GitLab to provide several features or customizations.

* [Description templates](https://docs.gitlab.com/ee/user/project/description_templates.html) for issues and merge requests
* define [code owners](https://docs.gitlab.com/ee/user/project/code_owners.html) for extended push/review rules
* for an extensive example see: https://gitlab.com/gitlab-org/gitlab/-/tree/master/.gitlab

We also collect some tools used in the CI pipeline here.
The [ci](ci) directory contains stuff used in the Windows builds.
