# AWS Changes

Account `217113048582` (iGEM Toronto), region `us-east-1`.
Access via AWS SSO profile `AdministratorAccess-217113048582`.

---

## 2026-08-02 — `esmfold2-centroids/60pid` 403s: not a permissions problem, the data was never uploaded

**No change applied — none is needed.** Reported as "public GET is broken on the `60pid` /
`60pid-msa` prefixes." It is not. The bucket policy already covers them, and the 403 is an artifact
of how S3 answers anonymous callers.

### The 403 means "no such key", not "access denied"

When a caller lacks `s3:ListBucket`, S3 deliberately returns **403 on a missing key** instead of
404 — otherwise anonymous probing would leak which keys exist. The public policy here grants
`GetObject` only (intentionally, see *Deliberately unchanged* below), so every nonexistent key is a
403 to the outside world.

Proof — a bogus key under the prefix that **does** work returns the same 403:

| Anonymous request | Result |
| --- | --- |
| `esmfold2-centroids/test2/structures/orf4981589.cif` | **200** |
| `esmfold2-centroids/test2/structures/THIS-KEY-DOES-NOT-EXIST.cif` | **403** |
| `esmfold2-centroids/60pid/structures/orf4981589.cif` | 403 |
| `esmfold2-centroids/60pid-msa/structures/orf4981589.cif` | 403 |

A working prefix and a broken one are indistinguishable from outside. The 403 carries no
information about permissions at all — chasing it as an ACL/BPA/encryption bug is a dead end.

### What is actually under those prefixes

Authenticated `HeadObject` as admin returns **404 Not Found** for both `60pid` and `60pid-msa` keys:

| Prefix | Objects | Size |
| --- | --- | --- |
| `esmfold2-centroids/60pid/` | 1 (zero-byte folder marker, created 2026-07-11) | 0 B |
| `esmfold2-centroids/90pid/` | 1 (zero-byte folder marker, created 2026-07-11) | 0 B |
| `esmfold2-centroids/60pid-msa/` | **does not exist** — not even a marker | — |
| `esmfold2-centroids/test/` | 21 | 4.0 MB |
| `esmfold2-centroids/test2/` | 21 | 5.0 MB |

`60pid/` and `90pid/` are empty placeholders someone made in the console three weeks before the
test runs. No `structures/` or `metrics/` subprefix exists under either. Nothing anywhere in the
bucket has `msa` in its key — `60pid-msa` has never existed under any spelling.

**Only `test/` and `test2/` (21 objects each) hold real folding output.** The real question is not
about S3; it is that the 60pid/60pid-msa ESMFold runs have not been uploaded.

### The hypotheses from the thread, all ruled out

| Suspected cause | Actual state |
| --- | --- |
| Block Public Access on the prefix | all four settings `false`; no account-level BPA |
| Object ACLs / cross-account ownership | `ObjectOwnership: BucketOwnerEnforced` — **ACLs are disabled bucket-wide**, so per-object ACLs cannot exist, let alone differ per prefix |
| Object encryption blocking anonymous reads | SSE-**AES256** (not KMS). The working `test2` object uses the identical setting; AES256 needs no key permission |
| Policy not covering the prefix | `Resource` is `arn:aws:s3:::petadex-protein-structures/*` — bucket-wide. There is no per-prefix grant to add |

The policy pasted in the thread had `"Principal": ""` and `Resource: .../` with the `*` characters
missing — a copy/paste artifact of Discord eating the asterisks. The live policy is correct.

**S3 prefixes are not directories and carry no permissions.** There is no such thing as opening
access "for `60pid/`" separately — anything written there is already public the moment it lands.

### Listing being 403 is correct and separate

Rya was right: listing is out of scope by design. `GetObject` is granted, `ListBucket` is not, so
the key list stays private. That is the deliberate shape (see below) and was not changed.

Consequence worth knowing: **there is no way for a collaborator to discover valid keys.** They must
be told exact keys, or be given a manifest. `test2/tracker.csv` is public and readable and appears
to serve that role:

```
https://petadex-protein-structures.s3.amazonaws.com/esmfold2-centroids/test2/tracker.csv
```

If the 60pid runs get uploaded, shipping a `tracker.csv` alongside them is the low-cost fix — it
keeps the bucket non-enumerable while making the contents discoverable.

### If you *do* want listing

Adding `s3:ListBucket` for `Principal: "*"` would make missing keys return an honest 404 and let
people browse. It also makes the whole bucket enumerable, which the current design deliberately
avoids and which raises the scraper/egress concern noted below. **Not applied — that is a call for
whoever owns the bucket, not a bug fix.** Scope it with a `s3:prefix` condition if you go ahead.

---

## 2026-08-02 — Research VMs had no credentials at all; gave them read-only S3

**Applied:** 14:56:37 EDT (18:56:37 UTC) by `igem-dennis`, associations immediately following
**Resources created:** IAM role `igem-vm-s3`, customer-managed policy `igem-vm-s3-access`,
instance profile `igem-vm-s3`; associated with **all 17 EC2 instances** in the account

### Problem

Reported as "the VMs can't write to S3". It was not a permissions problem — **the VMs had no
credential source whatsoever.**

| Check | Finding |
| --- | --- |
| Instance profiles attached | **0 of 17 instances** (`IamInstanceProfile: None`) |
| IAM users in `217113048582` | **zero** — the account cannot have issued static keys |
| Instances registered in SSM | **none** — corroborates the above; SSM registration itself requires an instance profile |
| IMDS config | enabled, IMDSv2 required, hop limit 2 — healthy, not the cause |
| Bucket policies containing `Deny` | none |
| SCPs | management account — never apply (same as the Bedrock entry below) |

With no role, IMDS serves no credentials and the SDK credential chain falls through entirely:
`NoCredentialsError` / "Unable to locate credentials".

Where VMs instead had `rnalab-*` keys from account `797308887321` baked into `~/.aws/credentials`,
that is a second, independent failure: those credentials get `AccessDenied`. Confirmed by direct
test against `petadex-signalp-300m-results`, `petadex-orf-fastaa`, and `petabite` — all denied. No
iGEM bucket grants cross-account access to `797308887321`, so those keys can never work here.

**Two red herrings, both ruled out:**

- The `DeleteUser` by root on 2026-06-08 was an **SSO Identity Store** user (directory
  `d-906670ea2a`, an offboarded teammate), *not* an IAM user. It did not break VM credentials.
- No `DisassociateIamInstanceProfile` events exist. Nothing was detached — the profiles were
  never there in the first place.

### Change

Scope was corrected mid-implementation: the role was first built read/write, then narrowed to
**read-only** before any instance was attached. **No VM ever held write permissions** — the
association step ran only after the narrowing. The write-capable policy version `v1` was
**deleted**, not merely superseded, so it cannot be reactivated by flipping the default version.

`igem-vm-s3-access` (v2, default and only version):

```json
{
  "Sid": "ReadOnlyAllOrgOwnedBuckets",
  "Effect": "Allow",
  "Action": ["s3:ListBucket","s3:GetBucketLocation","s3:GetObject","s3:GetObjectVersion"],
  "Resource": "*",
  "Condition": { "StringEquals": { "s3:ResourceAccount": "217113048582" } }
}
```

Plus `s3:ListAllMyBuckets` on `*` (unconditioned — the action has no resource to scope to).

**Why `ResourceAccount` instead of listing buckets:** covers buckets created later without a policy
edit, while still refusing to reach into other accounts. Enumerating the 17 current buckets would
silently fail to cover the 18th.

### Verification

`simulate-principal-policy` against `igem-vm-s3`:

| Simulated call | Decision |
| --- | --- |
| `s3:GetObject` on `petadex-orf-fastaa/*` | allowed |
| `s3:ListBucket` on `petadex-orf-fastaa` | allowed |
| `s3:PutObject` on `petadex-orf-fastaa/*` | implicitDeny |
| `s3:DeleteObject` on `petadex-orf-fastaa/*` | implicitDeny |

**Simulator gotcha — cost real debugging time.** `simulate-principal-policy` does **not**
auto-populate `s3:ResourceAccount`, so a correct policy reports `implicitDeny` on *reads* and looks
broken. The reads above only evaluate correctly with the key supplied explicitly:

```bash
--context-entries ContextKeyName=s3:ResourceAccount,ContextKeyValues=217113048582,ContextKeyType=string
```

Any conditioned policy needs this, or simulation results are meaningless. (Compare the note below:
simulation ignores SCPs, boundaries, and resource policies — this is a third blind spot.)

Association state: 17/17. The one running instance (`i-0724cb430ab5ff7b9`) reports `associated`;
the 16 stopped ones sit at `associating` until next boot. **That is normal and not a partial
failure** — the association is recorded and takes effect on start.

### This is a default, not a guardrail

**A signed-in user can still write.** Verified: the SSO `AdministratorAccess` role evaluates
`allowed` for `PutObject` and `DeleteObject`.

Two reasons, both load-bearing:

1. **Credential precedence.** The SDK checks env vars, `~/.aws/credentials`, and the SSO cache
   *before* IMDS. `aws sso login` on a VM puts the human ahead of the instance role, and the VM then
   acts as that human.
2. **Implicit vs explicit deny.** The read-only policy grants no writes, but that is an *implicit*
   deny — an absence of allow. It constrains only the principal it attaches to, and a signed-in user
   is a different principal entirely.

So this kills the accidental-overwrite class of problem for unattended jobs and anything using the
ambient instance identity. It is **not** an enforcement boundary against a human who signs in.

Enforcing that requires an **explicit `Deny` in a bucket policy** — the only lever that binds every
principal here, since SCPs don't apply to the management account. Not applied: an explicit Deny
would also break any current upload workflow, which is a decision for whoever owns those jobs.

### Operational caveat: stale keys shadow the role

If a VM still has `rnalab-*` keys in `~/.aws/credentials`, **attaching the role changes nothing** —
those keys keep winning by the same precedence rule and keep failing with `AccessDenied`. Clearing
them is required for reads to start working. Check with:

```bash
aws sts get-caller-identity     # expect .../igem-vm-s3, NOT rnalab-*
```

### Not applied

- **`AmazonSSMManagedInstanceCore` was not attached** — blocked by a local permission classifier.
  Not required for S3. Without it the VMs stay unreachable via Session Manager (SSH only), and they
  will continue to be absent from `describe-instance-information`.

### Left alone

`signalp-worker-role` — still exists, still attached to nothing. Even if attached it only permits
writes under `petadex-signalp-300m-results/signalp_other/{chunks,done,claims,failed}/*` (verified:
`signalp_other/done/x.txt` allowed, `other/x.txt` and `petadex-orf-fastaa/x.txt` both implicitDeny).
Not a general-purpose VM role; do not attach it expecting broad access.

---

## 2026-08-02 — Fix `bedrock-s3-write` InvokeModel grant

**Applied:** 14:45:44 EDT (18:45:44 UTC) by `igem-dennis`
**Resource:** IAM role `bedrock-s3-write`, inline policy `bedrock-s3-writePolicy`, statement `InvokeTheModel`

### Problem

The `bedrock:InvokeModel` grant named a model ARN that matched no model:

```
arn:aws:bedrock:us-east-1::foundation-model/anthropic.claude-4-6-sonnet-*
```

Two independent defects:

1. **Wrong ID ordering.** Bedrock IDs are `anthropic.` + the first-party model ID, which puts family
   before version: `anthropic.claude-sonnet-4-6`, not `claude-4-6-sonnet`. Verified against
   `aws bedrock list-foundation-models --by-provider anthropic`; nothing in the account's model list
   starts with `anthropic.claude-4-6-sonnet-`.

2. **Wrong resource type.** `anthropic.claude-sonnet-4-6` reports
   `inferenceTypesSupported: ["INFERENCE_PROFILE"]` with no `ON_DEMAND`. A bare foundation-model ARN
   is therefore not invocable on its own, so correcting the spelling alone would not have fixed it.

### Change

```json
{
  "Sid": "InvokeTheModel",
  "Effect": "Allow",
  "Action": ["bedrock:InvokeModel"],
  "Resource": [
    "arn:aws:bedrock:us-east-1:217113048582:inference-profile/us.anthropic.claude-sonnet-4-6",
    "arn:aws:bedrock:us-east-1::foundation-model/anthropic.claude-sonnet-4-6",
    "arn:aws:bedrock:us-east-2::foundation-model/anthropic.claude-sonnet-4-6",
    "arn:aws:bedrock:us-west-2::foundation-model/anthropic.claude-sonnet-4-6"
  ]
}
```

The two S3 statements (`ListTheBucket`, `ReadWriteObjects`) were preserved byte-identical.

### Why four ARNs

Invoking through an inference profile requires permission on **both** the profile ARN and the
underlying foundation-model ARNs in **every region the profile can route to**. The three
foundation-model entries are load-bearing, not redundant:

```
$ aws bedrock get-inference-profile --inference-profile-identifier us.anthropic.claude-sonnet-4-6
models: us-east-1, us-east-2, us-west-2
```

IAM is evaluated against whichever region actually served the request. A policy naming only the
profile plus `us-east-1` fails on a *fraction* of requests, which presents as flakiness rather than
misconfiguration.

**Scope decision:** chose the `us.` profile over `global.`. The `global.` variant routes to a
region-less ARN (`arn:aws:bedrock:::foundation-model/...`), i.e. inference may run outside the US.
`us.` keeps it in US jurisdiction, which is the right default for this dataset.

### Verification

`aws iam simulate-principal-policy` against the role:

| Simulated call | Decision |
| --- | --- |
| `bedrock:InvokeModel` on `us.anthropic.claude-sonnet-4-6` profile | allowed |
| `bedrock:InvokeModel` on old `claude-4-6-sonnet` style ARN | implicitDeny |
| `s3:PutObject` on `batch-extraction-.../out/...` | allowed (no regression) |

### Caller-side requirement

Pass the **profile ID**, not the bare model ID:

```python
modelId = "us.anthropic.claude-sonnet-4-6"    # correct
modelId = "anthropic.claude-sonnet-4-6"       # fails: requires an inference profile
```

The bare ID returns a *validation* error, not access-denied — so it will not look like a
permissions problem when it surfaces.

---

## 2026-08-02 — Batch job failure (resolved, no action taken)

Job `bio-extraction-1785695349`
(`arn:aws:bedrock:us-east-1:217113048582:model-invocation-job/0fzrtt7mb9xm`)
failed with `Customer doesn't have permissions to invokeModel`.

**Root cause: the job predates the policy fix by 8m22s.** Not a live misconfiguration.

| Time (EDT) | Event |
| --- | --- |
| 14:29:09 | Job submitted |
| 14:31–14:33 | Validating — passed (S3 read + JSONL format both OK) |
| 14:37:22 | Failed at end of `InProgress` |
| 14:45:44 | Policy fix applied |

The job's own record showed it was already calling the right things — `roleArn`
`bedrock-s3-write`, `modelId` the `us.anthropic.claude-sonnet-4-6` profile ARN. It ran against the
old broken resource. **Resubmitting unchanged is the fix.**

Failure point is diagnostic: failing during `Validating` is input-side (role can't read the input
prefix, malformed JSONL, or record count below the per-job minimum); failing at the end of
`InProgress` is model or output side.

Other causes ruled out:

- Model access — `authorizationStatus: AUTHORIZED`, `entitlementAvailability: AVAILABLE`
- No permissions boundary on the role
- SCPs — `217113048582` is the management account of org `o-u10flj6uvh`; SCPs never apply there
- Trust policy — correct (`bedrock.amazonaws.com` + `aws:SourceAccount` + `aws:SourceArn` scoped to
  `us-east-1` `model-invocation-job/*`)

Note that `simulate-principal-policy` evaluates the identity policy **only** — it does not consider
SCPs, permission boundaries, or resource policies. Passing simulation is necessary, not sufficient.

---

## Deliberately unchanged

### `petadex-protein-structures` is public — intended

```json
{"Sid":"PublicRead","Effect":"Allow","Principal":"*",
 "Action":"s3:GetObject","Resource":"arn:aws:s3:::petadex-protein-structures/*"}
```

`IsPublic: true`, all four Block Public Access settings `false`, no account-level BPA.

**This is the intended behaviour** — public structure downloads are the point of the database. Do
not "fix" this.

The policy shape is already the conservative one: `GetObject` only, no `ListBucket`, so the key list
stays private and the bucket is not enumerable.

Two consequences to keep in mind:

- **Anything written to this bucket becomes public by default.** Keep it to published structures only.
- **The real exposure is egress cost, not disclosure** — large public structure files plus a scraper
  is a billing story. Consider a CloudWatch billing alarm on `BytesDownloaded`. Requester Pays is the
  wrong tool here; it would break anonymous access.

Note: neither Bedrock role references this bucket, so it is unrelated to the batch pipeline.

### `enzyme-extraction-bedrock-batch-role` — known gaps, left as-is

- **No `bedrock:InvokeModel`.** Matches AWS's documented batch-inference service role, which is
  S3-only (the model permission is checked against the caller, not the service role). If this role is
  ever used for direct synchronous invocation, it needs the same `InvokeTheModel` block as
  `bedrock-s3-write`.
- **No `s3:GetBucketLocation`** (which `bedrock-s3-write` has). Cosmetic for batch.
- Uses `ArnEquals` where `bedrock-s3-write` uses `ArnLike`. **Not a defect** — AWS documents the two
  as behaving identically, both supporting `*` wildcards in ARN components.

---

## Open items

- [ ] **Upload the 60pid / 60pid-msa ESMFold output** — the reported "403 bug" is just missing data.
      `60pid/` and `90pid/` are empty folder markers; `60pid-msa/` does not exist. No S3 change needed.
- [ ] Decide whether to publish a `tracker.csv` manifest per prefix (as `test2/` has) so collaborators
      can discover keys without the bucket being made enumerable
- [ ] Resubmit `bio-extraction` batch job (unchanged input/role/modelId)
- [ ] Grep pipeline for a bare `anthropic.claude-sonnet-4-6` modelId; must be `us.anthropic.claude-sonnet-4-6`
- [ ] Add `bedrock:InvokeModelWithResponseStream` to `bedrock-s3-write` if anything streams
      (model reports `responseStreamingSupported: true`; only `InvokeModel` is currently granted)
- [ ] CloudWatch billing alarm on `BytesDownloaded` for `petadex-protein-structures`
- [ ] **No CloudTrail trail exists in this account at all** (`describe-trails` returns empty). No S3
      data events are recorded and only the 90-day Event History is queryable. Fix this *before*
      adding any Deny bucket policy — otherwise there is no way to tell what the Deny would break,
      and no audit trail of who wrote what. This also caps the forensics on the entry above.
- [ ] Clear stale `rnalab-*` keys from `~/.aws/credentials` on each VM; they shadow the new instance
      role. 16 of 17 are stopped, so this is per-VM work at next boot.
- [ ] Decide whether VM writes should be *enforced* off via explicit-Deny bucket policy, or whether
      "read-only by default, writes attributable to a signed-in human" is the desired posture
- [ ] Attach `AmazonSSMManagedInstanceCore` to `igem-vm-s3` (blocked locally this session) so the
      VMs become Session Manager reachable — currently SSH-only
- [ ] `igem-adi` called `PutRolePolicy` at 14:46:09 EDT, 25s after the fix. CloudTrail did not record
      `roleName`/`policyName` for that event. `bedrock-s3-write` was **not** clobbered (verified by
      direct read), but confirm what it touched — concurrent IAM edits will confound future debugging.

---

## Useful commands

```bash
export AWS_PROFILE=AdministratorAccess-217113048582

# Which Anthropic models exist here, and can they be invoked directly?
aws bedrock list-foundation-models --by-provider anthropic --region us-east-1 \
  --query 'modelSummaries[].{id:modelId,inference:inferenceTypesSupported}' --output table

# Which regions does a profile route to? (determines the IAM resource list)
aws bedrock get-inference-profile --inference-profile-identifier us.anthropic.claude-sonnet-4-6 \
  --region us-east-1 --query 'models[].modelArn'

# Why did a batch job fail? `message` carries the real reason.
aws bedrock get-model-invocation-job --region us-east-1 \
  --job-identifier <job-arn> \
  --query '{status:status,message:message,modelId:modelId,roleArn:roleArn}'

# Test an authorization decision without running a job (identity policy only)
aws iam simulate-principal-policy \
  --policy-source-arn arn:aws:iam::217113048582:role/bedrock-s3-write \
  --action-names bedrock:InvokeModel \
  --resource-arns arn:aws:bedrock:us-east-1:217113048582:inference-profile/us.anthropic.claude-sonnet-4-6

# Which VMs have an instance profile? (empty Profile column == no credentials at all)
aws ec2 describe-instances --query \
  'Reservations[].Instances[].{ID:InstanceId,Name:Tags[?Key==`Name`]|[0].Value,Profile:IamInstanceProfile.Arn}' \
  --output table

# Association state; stopped instances stay "associating" until next boot
aws ec2 describe-iam-instance-profile-associations \
  --query 'IamInstanceProfileAssociations[].{ID:InstanceId,State:State}' --output table

# Simulate a CONDITIONED policy — without --context-entries this reports a false implicitDeny
aws iam simulate-principal-policy \
  --policy-source-arn arn:aws:iam::217113048582:role/igem-vm-s3 \
  --action-names s3:GetObject s3:PutObject \
  --resource-arns arn:aws:s3:::petadex-orf-fastaa/x.txt \
  --context-entries ContextKeyName=s3:ResourceAccount,ContextKeyValues=217113048582,ContextKeyType=string \
  --query 'EvaluationResults[].{Action:EvalActionName,Decision:EvalDecision}' --output table

# Run ON the VM: confirms whether the instance role is actually in effect
aws sts get-caller-identity     # expect .../igem-vm-s3, NOT rnalab-*
```
