import { Moon } from "./icons";

export default function Nav({ onHome, onExample }: { onHome: () => void; onExample: (acc: string) => void }) {
  function toggleTheme() {
    const r = document.documentElement;
    const dark = r.getAttribute("data-theme") === "dark" ||
      (!r.getAttribute("data-theme") && matchMedia("(prefers-color-scheme: dark)").matches);
    r.setAttribute("data-theme", dark ? "light" : "dark");
  }
  return (
    <nav>
      <div className="wrap nav-in">
        <button className="brand" type="button" onClick={onHome} aria-label="Back to start">
          <span className="dot" />TmJunction
        </button>
        <div className="nav-links">
          <a href="#" className="keep-hide">How it works</a>
          <a href="#">Method</a>
          <a href="#">Docs</a>
          <span className="vdiv" />
          <button className="icon-btn" aria-label="Toggle theme" onClick={toggleTheme}><Moon /></button>
          <button className="btn" onClick={() => onExample("NM_001256799.3")}>Try GAPDH</button>
        </div>
      </div>
    </nav>
  );
}
