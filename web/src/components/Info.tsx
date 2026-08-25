import { useEffect, useRef, useState, type ReactNode } from "react";

/**
 * The one place explanatory prose lives now.
 *
 * Every panel used to carry its rationale inline, so the reading load landed on everyone
 * every time — including the people who already knew. The detail is still there, but it is
 * pulled rather than pushed: a short label states WHAT, this dot answers WHY on demand.
 *
 * Opens on hover for a quick look and on click/focus for a sticky one, so it is reachable
 * by keyboard and readable on touch, where hover does not exist.
 */
export default function Info({ children, label = "What is this?" }: {
  children: ReactNode;
  label?: string;
}) {
  const [open, setOpen] = useState(false);
  const [stuck, setStuck] = useState(false);   // clicked open: survives mouse-out
  const ref = useRef<HTMLSpanElement>(null);

  useEffect(() => {
    if (!stuck) return;
    const onKey = (e: KeyboardEvent) => { if (e.key === "Escape") { setStuck(false); setOpen(false); } };
    const onDown = (e: MouseEvent) => {
      if (!ref.current?.contains(e.target as Node)) { setStuck(false); setOpen(false); }
    };
    document.addEventListener("keydown", onKey);
    document.addEventListener("mousedown", onDown);
    return () => {
      document.removeEventListener("keydown", onKey);
      document.removeEventListener("mousedown", onDown);
    };
  }, [stuck]);

  return (
    <span className="info" ref={ref}
      onMouseEnter={() => setOpen(true)}
      onMouseLeave={() => { if (!stuck) setOpen(false); }}>
      <button type="button" className="info-dot" aria-label={label} aria-expanded={open}
        onClick={(e) => { e.stopPropagation(); setStuck((v) => !v); setOpen((v) => !v || !stuck); }}
        onFocus={() => setOpen(true)}
        onBlur={() => { if (!stuck) setOpen(false); }}>i</button>
      {open && <span className="info-pop" role="tooltip">{children}</span>}
    </span>
  );
}
