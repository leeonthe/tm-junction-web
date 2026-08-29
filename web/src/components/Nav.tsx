import { useEffect, useState } from "react";
import { Moon, Sun } from "./icons";

const THEME_KEY = "tmj.theme";
type Theme = "light" | "dark";

const systemTheme = (): Theme =>
  matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light";

/** What is on screen right now: an explicit choice if one was made, else the OS setting. */
function currentTheme(): Theme {
  const set = document.documentElement.getAttribute("data-theme");
  return set === "dark" || set === "light" ? set : systemTheme();
}

export default function Nav({ onHome, onMethod, methodOn }: {
  onHome: () => void; onMethod: () => void; methodOn: boolean;
}) {
  const [theme, setTheme] = useState<Theme>(currentTheme);

  // A toggle that forgets on reload reads as unfinished. The choice is stored, and
  // index.html applies it before first paint so the page never flashes the wrong theme.
  function toggleTheme() {
    const next: Theme = theme === "dark" ? "light" : "dark";
    document.documentElement.setAttribute("data-theme", next);
    try { localStorage.setItem(THEME_KEY, next); } catch { /* private mode */ }
    setTheme(next);
  }

  // Follow the OS while the user has not overridden it.
  useEffect(() => {
    const mq = matchMedia("(prefers-color-scheme: dark)");
    const sync = () => { if (!localStorage.getItem(THEME_KEY)) setTheme(systemTheme()); };
    mq.addEventListener("change", sync);
    return () => mq.removeEventListener("change", sync);
  }, []);

  return (
    <nav aria-label="Main">
      <div className="wrap nav-in">
        <button className="brand" type="button" onClick={onHome}
          aria-label="Exon Junction Primer — start a new search">
          <Mark />
          <span className="brand-name">Exon Junction Primer</span>
        </button>

        <div className="nav-links">
          <button type="button" className={`navlink ${methodOn ? "on" : ""}`}
            aria-current={methodOn ? "page" : undefined} onClick={onMethod}>Method</button>
          <button className="icon-btn" type="button" onClick={toggleTheme}
            aria-label={`Switch to ${theme === "dark" ? "light" : "dark"} theme`}>
            {theme === "dark" ? <Sun /> : <Moon />}
          </button>
        </div>
      </div>
    </nav>
  );
}

/** Two exons joined over a spliced-out intron — the diagram the whole tool is about. */
function Mark() {
  return (
    <svg className="brand-mark" viewBox="0 0 24 24" aria-hidden="true">
      <rect width="24" height="24" rx="7" fill="var(--brand)" />
      <path d="M10.4 13.4 12 8.6l1.6 4.8" fill="none" stroke="#fff" strokeWidth="1.7"
        strokeLinecap="round" strokeLinejoin="round" />
      <rect x="3.6" y="13.2" width="6.8" height="5.6" rx="1.6" fill="#fff" />
      <rect x="13.6" y="13.2" width="6.8" height="5.6" rx="1.6" fill="#fff" />
    </svg>
  );
}
