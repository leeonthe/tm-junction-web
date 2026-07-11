/**
 * DNA double-helix loader — the original DNA helix art, reduced to 6 rungs, recolored to
 * brand blue #237AF2 (near) / grey (far), hue-cycling removed, `d`-morph motion preserved.
 * Embedded via <object> so the SMIL animation plays. `horizontal` uses the 90°-rotated art.
 */
const V_ASPECT = 0.705;   // dna-blue.svg viewBox w/h (portrait)

export default function DnaLoader({
  size = 150, horizontal = false, label,
}: {
  size?: number;        // the box's long side in px
  horizontal?: boolean;
  label?: string;
}) {
  const w = horizontal ? size : Math.round(size * V_ASPECT);
  const h = horizontal ? Math.round(size * V_ASPECT) : size;
  const src = horizontal ? "/dna-blue-h.svg" : "/dna-blue.svg";
  return (
    <div className="dna-loader">
      <object className="dna-svg" type="image/svg+xml" data={src}
        style={{ width: w, height: h }} aria-label={label || "Loading"} tabIndex={-1}>
        DNA
      </object>
      {label && <div className="dna-label">{label}</div>}
    </div>
  );
}
