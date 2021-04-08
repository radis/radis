# Contributing to RADIS

We love spectroscopy, and we love your input! We want to make contributing to this project as easy as possible, whether it's:

- Reporting a bug
- Discussing the current state of the code or the physics
- Submitting a fix
- Proposing new features (interface or physics)
- Becoming a core developer

This document is a quick summary, more information can be found on the [Developer Guide](https://radis.readthedocs.io/en/latest/dev/developer.html)


## Report bugs using Github's [issues](https://github.com/radis/radis/issues)
We use GitHub issues to track public bugs. Report a bug by [opening a new issue](https://github.com/radis/radis/issues/new/choose); it's that easy!

**Great Bug Reports** tend to have:

- A quick summary and/or background
- Steps to reproduce : specific, with code samples if possible
- What you expected would happen
- What actually happens
- Notes (possibly including why you think this might be happening, or stuff you tried that didn't work)

Issues labeled [Good First Issue](https://github.com/radis/radis/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22)
are also a nice way to get started with contributing to the project.


## Ask questions

Join the community chat:

[![Slack](https://img.shields.io/badge/slack-join-green.svg?logo=slack)](https://radis.github.io/slack-invite/)
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/radis-radiation/community)

## Propose new features

You can suggest or vote for new features below:

[![Feature Requests](https://feathub.com/radis/radis?format=svg)](https://feathub.com/radis/radis)


## Change the code : become a Contributor.

We use the [Github Flow](https://guides.github.com/introduction/flow/index.html), where
all code changes should happen through Pull Requests (PR). PR are the best way to propose changes to the codebase :
Maintainers of the code can comment, help you, and even edit your code directly.

1. **Fork** the repo : this creates your own RADIS version.
2. If it's an easy fix (ex : fix a typo, update a single file, ...), edit it online on your Fork then open the PR (go directly to 6.)
3. Else, **Clone** your fork : this adds RADIS on your local computer.

```git clone https://github.com/[YOUR-USERNAME]/radis.git```

4. Time to work now : make your changes locally ! If you plan to work on multiple fixes at the same time,
create a new branch (with a descriptive name) for your feature.

5. Push your changes and open a Pull request (PR). If things aren't final, it's okay: just mark is as a draft `[WIP]`. Maintainers will review and start helping you from there !
6. Once the review is complete (physical tests + code linting), the pull request is merged : welcome to the [RADIS contributors](https://github.com/radis/radis/graphs/contributors) ! :clap:

*TODO : add a .gif for the SmartGit/GitHub Desktop version ?*

## Where to start ?

If it's your first time contributing to an open-source project, welcome ! 👋

The best is to get familiar with the procedure above. For instance, there are many tiny improvements to be made to the Documentation :
- Have a look at the [Documentation TODO List](https://github.com/radis/radis/issues/77).
- Pick one of them, create your GitHub account, and start the procedure above to become a Contributor.

Then, have a look at the [GitHub opened issues](https://github.com/radis/radis/issues). Easy issues that will help you understand the code structure are labelled as
[Good First Issue](https://github.com/radis/radis/contribute).


## Regular Contributor ? Keep your Fork updated.

If you keep on using your own local RADIS version, you want to keep it updated with the main branch.
In usual Git Flow, your own Fork ``[YOUR-USERNAME]/radis`` is refered to as ``origin``, and the main
repo ``radis/radis`` as ``upstream``.

One workflow is to have at least 2 branches locally :

- one pointing to the latest developer version, e.g. named, ``upstream-develop``.
- one per feature or fix you're currently working on, poiting to your own fork ``origin`` : could be ``develop``, ``fix/something``, ``this_new_idea_i_work_on`` , etc.

A way to do this is to :

1. set-up another remote :

```
git remote add upstream git://github.com/radis/radis.git
git fetch upstream
```

2. Starting a new fix/feature ? Branch from the latest ``upstream-develop``.
```
git branch -b [NEW_BRANCH] upstream/develop
```

3. Need to update your local branch, by 'rebasing', i.e. without creating a merge-commit for the changes that were done by others  :

```
git pull --rebase upstream [NEW_BRANCH]
```
4. Push your changes to your own Fork
```
git push -u origin [NEW-BRANCH]
```
5. Open a PR, as in "Become a Contributor" above.



## Linting : use a Consistent Coding Style

We're using the Black coding style, and Flake8 to ensure no syntax errors are introduced.
If you're a first time contributor, do not worry about Linting. Other developers will do it for you.

Once you're a regular contributor, we advise you install pre-commit that will take care of all the formatting :

```
cd radis
pre-commit install
```

See the [Developer Guide](https://radis.readthedocs.io/en/latest/dev/developer.html#code-linting)
for more information.

## License
By contributing, you agree that your contributions will be licensed under the [GNU Lesser General Public License v3.0](https://github.com/radis/radis/blob/develop/LICENSE)

## References
This document was adapted from the open-source contribution guidelines of [@briandk](https://gist.github.com/briandk/3d2e8b3ec8daf5a27a62)
