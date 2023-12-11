# Notes for developers

## Tagging a new version

This reusable action follows [Semantic Versioning, v2.0.0][1].

A new __patch__ version is automatically tagged when a group of commits is pushed to the `main` branch;
for example, when a group that comprises a pull request is merged.
Alternatively, a new patch version is tagged for each commit in the group that has a message title prefixed with `fix`.
For example, a commit with the following message title would tag a new patch version when it is pushed to the `main` branch:

```
fix: a bug fix
```

A new __minor__ version is tagged for each commit in the group that has a message title prefixed with `feat`.
For example, a commit with the following message title would tag a new minor version when it is pushed to the `main` branch:

```
feat: a new feature
```

A new __major__ version is tagged for each commit in the group that has `BREAKING CHANGE` in its message body.
For example, a commit with the following message body would tag a new major version:

```
Remove a function

BREAKING CHANGE: Removing a function is not backwards-compatible.
```

Whilst there are other prefixes besides `fix` and `feat`, they do not tag new versions.

[1]: https://semver.org/spec/v2.0.0.html
