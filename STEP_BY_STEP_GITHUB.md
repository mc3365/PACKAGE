# Setting up `PACKAGE` on your laptop and GitHub — quick guide

This guide gets the skeleton onto your laptop, installs it, and pushes it to a new
private GitHub repo. You'll do most of this through **Claude Code** running locally, so
this guide is short — the heavy lifting is interactive.

Estimated time: **30-60 minutes**, depending on how much is already installed.

---

## Step 0 — One-time prerequisites

Before starting, make sure you have:

- [ ] **Miniconda installed** on your laptop. Download from
      https://docs.conda.io/en/latest/miniconda.html and follow the installer. ~5 min.
- [ ] **Git installed.** Check with `git --version`. macOS usually has it via Xcode
      Command Line Tools. Linux: `apt install git` or similar.
- [ ] **Claude Code installed** on your laptop. See https://claude.com/claude-code
      for install instructions. Run `claude --version` to verify.
- [ ] **GitHub account.** You said you have this.

You do **not** need to set up GitHub SSH keys ahead of time — Claude Code will help
with that when needed.

---

## Step 1 — Unpack the skeleton

The skeleton folder I gave you is currently named `PACKAGE`. Move it to wherever you
want to keep code projects on your laptop (e.g. `~/projects/`):

```bash
# Adjust paths as needed
mv ~/Downloads/PACKAGE /Users/cmy324/Desktop/analysis/git/fiberforge/PACKAGE
cd /Users/cmy324/Desktop/analysis/git/fiberforge/PACKAGE
```

---

## Step 2 — Pick a name (optional, can defer)

Throughout the code, the placeholder name `PACKAGE` is used. You can keep this for now
and rename later when you've decided. To rename:

```bash
# Replace PACKAGE everywhere
find . -type f \( -name "*.py" -o -name "*.toml" -o -name "*.md" -o -name "*.yml" \
    -o -name "*.yaml" -o -name "*.cff" \) -not -path "./.venv/*" -exec \
    sed -i.bak 's/PACKAGE/your_chosen_name/g' {} +

# Rename the package directory
mv src/PACKAGE src/your_chosen_name

# Clean up backup files sed created
find . -name "*.bak" -delete
```

I recommend keeping `PACKAGE` for now and renaming once you're sure of the name.

---

## Step 3 — Hand off to Claude Code

Open a terminal in the `PACKAGE` folder and start Claude Code:

```bash
cd ~/projects/PACKAGE
claude
```

Then paste this initial prompt to give Claude Code the context:

> I have a Python package skeleton in this directory generated from a previous chat.
> I need help setting it up locally and pushing to a new private GitHub repo. Please:
>
> 1. Create a fresh conda environment with Python 3.11 named `PACKAGE_env`
> 2. Activate it and install the package with `pip install -e ".[dev]"`
> 3. Verify the install by running `PACKAGE --version` and `pytest`
> 4. If anything fails, debug it
> 5. Once tests pass, help me create a new **private** GitHub repository named
>    `PACKAGE` and push this code to it
> 6. Verify GitHub Actions CI runs and passes on the first commit
>
> The skeleton follows the standard src-layout. Tests should all pass already (they
> test placeholder functionality). After this step the package will have stub modules
> ready to be filled in.

Claude Code will walk you through it step by step. It will:

- Run installs and report progress
- Show you errors if anything breaks, then propose fixes
- Walk you through GitHub auth setup if needed (SSH key or personal access token)
- Help create the repo and push
- Confirm CI ran successfully

---

## Step 4 — Verify CI passed

After pushing, open the repo on github.com in your browser:

1. Go to **Actions** tab
2. You should see one workflow run for the initial commit
3. It will run a matrix of 2 jobs (Python 3.10 and 3.12 on Ubuntu)
4. Both should complete with a green checkmark within ~3 minutes

If anything goes red, copy the error log and ask Claude Code to fix it.

---

## Step 5 — You're done with Phase 1

You now have:

- ✅ A working local install on your laptop
- ✅ A private GitHub repository with the skeleton
- ✅ CI passing on every push
- ✅ A scaffold ready to fill in module by module

Come back to the main conversation (or continue with Claude Code) when you're ready
for **Phase 2: porting the database query API**. That's where we start filling in
real code, starting from your existing `fiber_database_v8.py`.

---

## Troubleshooting (if Claude Code can't fix it interactively)

**`claude` command not found**

You haven't installed Claude Code yet. See https://claude.com/claude-code.

**`conda` command not found after install**

You may need to close and reopen your terminal, or run `source ~/.bashrc`
(or `~/.zshrc` on macOS).

**Permission denied (publickey) when pushing to GitHub**

Either set up SSH keys (https://docs.github.com/en/authentication/connecting-to-github-with-ssh)
or use the HTTPS URL with a personal access token
(https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens).
Claude Code can walk you through either path.

**GitHub Actions failing**

If CI is red, click into the failing job on GitHub, copy the error log, and paste it
into Claude Code. Most likely cause for first-time CI failures is a Python version or
dependency issue that we can fix in `pyproject.toml`.
