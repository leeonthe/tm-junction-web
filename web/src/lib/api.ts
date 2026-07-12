import type { AnalyzeResponse, ApiError } from "./types";

const API_URL = (import.meta.env.VITE_API_URL as string | undefined) ?? "http://localhost:8000";

export class AnalyzeError extends Error {
  code: string;
  constructor(code: string, message: string) {
    super(message);
    this.code = code;
  }
}

export interface RemoteSuggestion { accession: string; gene: string }

export async function suggest(q: string, signal?: AbortSignal): Promise<RemoteSuggestion[]> {
  try {
    const res = await fetch(`${API_URL}/suggest?q=${encodeURIComponent(q)}`, { signal });
    if (!res.ok) return [];
    const d = await res.json();
    return (d.suggestions ?? []) as RemoteSuggestion[];
  } catch {
    return [];   // network/abort — fall back to local matches only
  }
}

export interface Progress { pct: number; detail: string }

/**
 * Streaming analysis via SSE — reports REAL progress (mostly the per-isoform sequence
 * fetches) and resolves with the result. Parses the event stream from a fetch body reader.
 */
export async function analyzeStream(
  accession: string,
  onProgress: (p: Progress) => void,
  signal?: AbortSignal,
): Promise<AnalyzeResponse> {
  let res: Response;
  try {
    res = await fetch(`${API_URL}/analyze/stream?accession=${encodeURIComponent(accession)}`, { signal });
  } catch {
    throw new AnalyzeError("NETWORK", "Can't reach the engine. Is the Python backend running?");
  }
  if (!res.ok || !res.body) {
    throw new AnalyzeError("NETWORK", "Can't reach the engine. Is the Python backend running?");
  }
  const reader = res.body.getReader();
  const decoder = new TextDecoder();
  let buf = "";
  for (;;) {
    const { done, value } = await reader.read();
    if (done) break;
    buf += decoder.decode(value, { stream: true });
    let sep: number;
    while ((sep = buf.indexOf("\n\n")) >= 0) {
      const frame = buf.slice(0, sep);
      buf = buf.slice(sep + 2);
      let event = "message";
      let data = "";
      for (const line of frame.split("\n")) {
        if (line.startsWith("event:")) event = line.slice(6).trim();
        else if (line.startsWith("data:")) data += line.slice(5).trim();
      }
      if (!data) continue;
      const parsed = JSON.parse(data);
      if (event === "progress") onProgress(parsed as Progress);
      else if (event === "result") return parsed as AnalyzeResponse;
      else if (event === "failed") throw new AnalyzeError(parsed.error, parsed.message);
    }
  }
  throw new AnalyzeError("ERROR", "The analysis stream ended unexpectedly.");
}

export async function analyze(accession: string): Promise<AnalyzeResponse> {
  let res: Response;
  try {
    res = await fetch(`${API_URL}/analyze`, {
      method: "POST",
      headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ accession }),
    });
  } catch {
    throw new AnalyzeError("NETWORK", "Can't reach the engine. Is the Python backend running?");
  }
  const data = (await res.json()) as AnalyzeResponse | ApiError;
  if (!res.ok || "error" in data) {
    const err = data as ApiError;
    throw new AnalyzeError(err.error ?? "ERROR", err.message ?? "Analysis failed.");
  }
  return data as AnalyzeResponse;
}
