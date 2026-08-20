import reactHooks from "eslint-plugin-react-hooks";
import tseslint from "typescript-eslint";

/**
 * Deliberately narrow: the ONE rule that catches a class of bug TypeScript cannot.
 *
 * A hook called after an early `return` changes the hook count between renders, and React
 * responds by throwing "Rendered more hooks than during the previous render" — which
 * unmounts the whole tree and leaves a blank page. That is not a subtle degradation the
 * user can work around; it is the app disappearing. It happened for real: JunctionDesigner
 * returned early for a transcript with no EEJ design, so switching from such a transcript
 * to one WITH a design (routine on a multi-isoform gene like TCF7L2) added a hook and blew
 * up the page, and only a reload recovered it.
 *
 * Style and exhaustive-deps rules are intentionally NOT enabled: this config exists to stop
 * the app from crashing, and a wall of pre-existing warnings would bury exactly that signal.
 */
export default tseslint.config(
  { ignores: ["dist/**", "node_modules/**"] },
  {
    files: ["src/**/*.{ts,tsx}"],
    // exhaustive-deps is off here, so its existing disable comments are not "unused" —
    // they document intent for whoever turns the rule on. Don't report them as problems.
    linterOptions: { reportUnusedDisableDirectives: "off" },
    languageOptions: { parser: tseslint.parser },
    plugins: { "react-hooks": reactHooks },
    rules: { "react-hooks/rules-of-hooks": "error" },
  },
);
