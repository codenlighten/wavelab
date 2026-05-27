# wavelab — deploy keys

Public SSH keys on the GPU droplet (`159.89.127.151` / `Ramp`,
RTX 6000 Ada). Add the one you want to give wavelab access on
GitHub → `codenlighten/wavelab` → **Settings → Deploy keys → Add**.

## Recommended for wavelab

The `rtx-server-tests-...` key is the most clearly scoped name for
this work. Give it **Allow write access** if you want droplet-side
commits to push back to the repo; otherwise leave write off for
pull-only mirroring.

```
ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIAfVvt4OYupkwnHmwaQ6lRBIlzLtw2+CKzNKaISWdWI6 rtx-server-tests-20260520
```

## All keys present on the droplet

| Label | Key |
| --- | --- |
| `metanation-droplet` | `ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIFV6cUBgbbtmxioPnpG8CdZeo5cIyq3qzALM9/R30NVa metanation-droplet` |
| `dsc2-deploy` | `ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIM3Our2IUrsTg/7bvjYbV9e988TdkLB/WtXns93XXOBu dsc2-deploy` |
| `rtx-server-tests-20260520` | `ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIAfVvt4OYupkwnHmwaQ6lRBIlzLtw2+CKzNKaISWdWI6 rtx-server-tests-20260520` |

Each of these is already in `/root/.ssh/authorized_keys`-style usage
on the droplet (they're listed in `~/.ssh/id_*.pub`). One of them is
also wired into GitHub as `OriginNeuralAI` — that account does NOT
have access to `codenlighten/wavelab`, which is why the initial
`git clone` attempt failed.

## After you add a deploy key

On the droplet:

```sh
cd /root/wavelab/wavelab
git remote -v                          # should already show origin
git fetch origin                        # confirms read access
# git push origin main                  # only if you gave the key write
```

If you'd rather not use deploy keys, two other paths work:

- **Make the repo public** — drop the auth requirement entirely.
- **Personal access token** — `git remote set-url origin https://<token>@github.com/codenlighten/wavelab.git`.
