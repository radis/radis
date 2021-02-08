# Contributing to RADIS

We love spectroscopy, and we love your input! We want to make contributing to this project as easy as possible, whether it's:

- Reporting a bug
- Discussing the current state of the code or the physics
- Submitting a fix
- Proposing new features (interface or physics)
- Becoming a core developer

This document is a quick summary, more information can be found on the [Developer Guide](https://radis.readthedocs.io/en/latest/dev/developer.html)

## We Develop with Github
We use github to host code, to track issues and feature requests, as well as accept pull requests.

## We (try to) use [Github Flow](https://guides.github.com/introduction/flow/index.html)
So all code changes (should) happen through Pull Requests.
Pull requests are the best way to propose changes to the codebase:

1. Fork the repo and create your branch from `develop`.
2. Clone your fork, make your changes locally (for a quick change you may edit it directly online with the GitHub editor).
3. Push your changes and open a Pull request (PR). If things aren't final, it's okay: just mark is as a draft `[WIP]`.
4. PR are a good support for discussion on your changes. Maintainers of the code can comment, help you, and even edit your code directly.
5. If you've added code that should be tested, add tests.
6. Once the physical tests pass, lint your code.
7. The pull request is merged, welcome to the [RADIS contributors](https://github.com/radis/radis/graphs/contributors)

## Report bugs using Github's [issues](https://github.com/radis/radis/issues)
We use GitHub issues to track public bugs. Report a bug by [opening a new issue](); it's that easy!

Issues labeled [Good First Issue](https://github.com/radis/radis/issues?q=is%3Aissue+is%3Aopen+label%3A%22good+first+issue%22)
are also a nice way to get started with contributing to the project.

## Write bug reports with detail, background, and sample code

**Great Bug Reports** tend to have:

- A quick summary and/or background
- Steps to reproduce
  - Be specific!
  - Give sample code if you can.
- What you expected would happen
- What actually happens
- Notes (possibly including why you think this might be happening, or stuff you tried that didn't work)

People *love* thorough bug reports. Not even kidding.

## Use a Consistent Coding Style

We're using the Black coding style. See the [Developer Guide](https://radis.readthedocs.io/en/latest/dev/developer.html#code-linting)
for more information. Can you also install pre-commit that will take care of all the formatting :

```
cd radis
pre-commit install
```

## Ask questions

Join the community chat:

[![Slack](https://img.shields.io/badge/slack-join-green.svg?logo=slack)](https://radis.github.io/slack-invite/)
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/radis-radiation/community)

## Propose new features

You can suggest or vote for new features below:

[![Feature Requests](https://feathub.com/radis/radis?format=svg)](https://feathub.com/radis/radis)

## License
By contributing, you agree that your contributions will be licensed under the [GNU Lesser General Public License v3.0](https://github.com/radis/radis/blob/develop/LICENSE)

## References
This document was adapted from the open-source contribution guidelines of [@briandk](https://gist.github.com/briandk/3d2e8b3ec8daf5a27a62)
