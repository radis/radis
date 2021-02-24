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

If it's your first time contributing to an open-source project, have a look at [Github Flow](https://guides.github.com/introduction/flow/index.html)

All code changes should happen through Pull Requests (PR). PR are the best way to propose changes to the codebase :
Maintainers of the code can comment, help you, and even edit your code directly.

1. **Fork** the repo : this creates your own RADIS version.
2. If it's an easy fix (ex : fix a typo, update a single file, ...), edit it online on your Fork then open the PR (go directly to 6.)
3. Else, **Clone** your fork : this adds RADIS on your local computer.

```git clone https://github.com/[MY-NAME]/radis.git```

4. Time to work now : make your changes locally ! If you plan to work on multiple fixes at the same time,
create a new branch (with a descriptive name) for your feature.

5. Push your changes and open a Pull request (PR). If things aren't final, it's okay: just mark is as a draft `[WIP]`. Maintainers will review and start helping you from there !
6. Once the review is complete (physical tests + code linting), the pull request is merged : welcome to the [RADIS contributors](https://github.com/radis/radis/graphs/contributors) ! :clap:


# TODO : add a gif for the SmartGit version ?


## Change the code : keep your Fork updated.

If you keep on using your own local RADIS version, you want to keep it updated with the main branch.
- The idea is

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
