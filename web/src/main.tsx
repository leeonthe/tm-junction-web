import React from "react";
import ReactDOM from "react-dom/client";
import App from "./App";
import TestPage from "./TestPage";
import "./styles/app.css";
import { trackVisit } from "./lib/track";

// Minimal routing: /test previews components; everything else is the app.
const isTest = window.location.pathname.replace(/\/+$/, "") === "/test";

// One visit per page load; the server keeps one per visitor per day (see lib/track).
if (!isTest) trackVisit();

ReactDOM.createRoot(document.getElementById("root")!).render(
  <React.StrictMode>
    {isTest ? <TestPage /> : <App />}
  </React.StrictMode>
);
