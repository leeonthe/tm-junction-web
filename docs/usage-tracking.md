# Usage tracking

TransPrimer records two things: who visited each day, and which genes or accessions were
searched. Code: [`web/src/lib/track.ts`](../web/src/lib/track.ts) (browser) and
[`web/api/_track.py`](../web/api/_track.py) (the `/api/track` endpoint).

## Turning it on

Set these on the Vercel project (web), then redeploy:

| Variable | Value |
|---|---|
| `MONGODB_URI` | the MongoDB connection string |
| `MONGODB_DB` | optional database name; defaults to `transprimer` |

Without `MONGODB_URI` the endpoint accepts events and stores nothing. Local dev servers
never send events; only production builds do. Automated browsers and crawlers are skipped.

## What is stored

A visitor is a random id kept in the browser's `localStorage` (`tmj.visitor`). It marks a
browser, not a person. `day` is the date in Eastern time (`America/New_York`); `at` and
`first_at` are UTC.

**`visits`**: one document per visitor per day.

```json
{ "visitor": "cd2b89e9-…", "day": "2026-11-19", "first_at": "2026-11-19T16:06:16Z" }
```

**`searches`**: one document per search the user started.

```json
{ "visitor": "cd2b89e9-…", "day": "2026-11-19", "at": "2026-11-19T16:06:24Z",
  "kind": "accession", "query": "NM_002046.7", "species": "human", "gene": "GAPDH", "ok": true }
```

A failed search has no `gene`, and adds what is needed to re-run it:

```json
{ "visitor": "cd2b89e9-…", "day": "2026-11-19", "at": "2026-11-19T16:06:45Z",
  "kind": "gene", "query": "Gapdh1x", "species": "fly", "ok": false,
  "error": "NOT_FOUND", "message": "No fruit fly (Drosophila melanogaster) gene found for “Gapdh1x”.",
  "engine": "same-origin" }
```

- `query` is what the user submitted, cut at 200 characters.
- `engine` is `configured` (the `VITE_API_URL` engine) or `same-origin` (the copy under `/api`).
- A failed accession search has no `species`: it never resolved to one.

**Counted as a search:** a gene search, an accession search, picking a transcript from the
gene's variant list, opening a shared link, and an "Existing" row in Custom sequence mode.

**Not counted:** back/forward presses, switching to a sibling transcript inside a result,
changing the excluded transcripts, and names the Custom row looks up while the user is
still typing (unless they resolve).

## Queries

Run in the Atlas Data Explorer aggregation tab or `mongosh`.

Visitors per day:

```js
db.visits.aggregate([{ $group: { _id: "$day", visitors: { $sum: 1 } } }, { $sort: { _id: -1 } }])
```

Visitors on one day:

```js
db.visits.countDocuments({ day: "2026-11-19" })
```

Returning visitors (seen on more than one day):

```js
db.visits.aggregate([
  { $group: { _id: "$visitor", days: { $sum: 1 }, first: { $min: "$day" }, last: { $max: "$day" } } },
  { $match: { days: { $gt: 1 } } },
  { $sort: { days: -1 } },
])
```

Top genes, with how many different visitors searched each:

```js
db.searches.aggregate([
  { $match: { ok: true } },
  { $group: { _id: { gene: "$gene", species: "$species" }, searches: { $sum: 1 }, visitors: { $addToSet: "$visitor" } } },
  { $project: { searches: 1, visitors: { $size: "$visitors" } } },
  { $sort: { visitors: -1, searches: -1 } },
  { $limit: 50 },
])
```

Searches by species:

```js
db.searches.aggregate([{ $group: { _id: "$species", searches: { $sum: 1 } } }, { $sort: { searches: -1 } }])
```

Failed searches, most frequent first:

```js
db.searches.aggregate([
  { $match: { ok: false } },
  { $group: { _id: { query: "$query", species: "$species", error: "$error" }, times: { $sum: 1 }, last: { $max: "$at" }, engines: { $addToSet: "$engine" } } },
  { $sort: { times: -1 } },
])
```

Everything one visitor searched:

```js
db.searches.find({ visitor: "cd2b89e9-…" }).sort({ at: 1 })
```

## Tests

```sh
engine/.venv/bin/python -m pytest web/api/_tests   # endpoint
cd web && npm test                                 # browser side (lib/track.test.ts)
```
