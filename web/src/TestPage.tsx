import DnaLoader from "./components/DnaLoader";

/** Preview page for the loading animation — reachable at /test. */
export default function TestPage() {
  return (
    <div className="test-page">
      <div className="test-head">
        <h1>Loading animation</h1>
        <a className="btn" href="/">Back to app</a>
      </div>
      <p className="test-sub">Original DNA helix — 6 rungs, brand blue, no hue cycling. Motion preserved. Vertical &amp; horizontal.</p>

      <div className="test-grid">
        <div className="test-cell">
          <DnaLoader size={180} label="Analyzing…" />
          <span className="test-cap">Vertical · large</span>
        </div>
        <div className="test-cell">
          <DnaLoader size={110} />
          <span className="test-cap">Vertical · small</span>
        </div>
        <div className="test-cell">
          <DnaLoader size={240} horizontal label="Analyzing…" />
          <span className="test-cap">Horizontal · large</span>
        </div>
        <div className="test-cell">
          <DnaLoader size={150} horizontal />
          <span className="test-cap">Horizontal · small</span>
        </div>
      </div>

      <p className="test-note">Usage: <code>{`<DnaLoader size={140} horizontal label="Analyzing…" />`}</code></p>
    </div>
  );
}
