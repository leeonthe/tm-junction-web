import React from "react";
import ReactDOM from "react-dom/client";
import App from "./App";
import TestPage from "./TestPage";
import "./styles/app.css";

// Minimal routing: /test previews components; everything else is the app.
const isTest = window.location.pathname.replace(/\/+$/, "") === "/test";

ReactDOM.createRoot(document.getElementById("root")!).render(
  <React.StrictMode>
    {isTest ? <TestPage /> : <App />}
  </React.StrictMode>
);
