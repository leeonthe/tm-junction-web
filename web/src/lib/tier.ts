import type { Tier } from "./types";

export const tierColorVar: Record<Tier, string> = {
  CONVENTIONAL: "var(--conv)",
  NEEDS_EEJ: "var(--eej)",
  NO_SINGLE_UNIQUE_JUNCTION: "var(--hard)",
};

export const tierChipClass: Record<Tier, string> = {
  CONVENTIONAL: "tc-conv",
  NEEDS_EEJ: "tc-eej",
  NO_SINGLE_UNIQUE_JUNCTION: "tc-hard",
};

export const tierLabel: Record<Tier, string> = {
  CONVENTIONAL: "Conventional",
  NEEDS_EEJ: "Needs EEJ",
  NO_SINGLE_UNIQUE_JUNCTION: "Hard case",
};

export const verdictClass: Record<Tier, string> = {
  CONVENTIONAL: "t-conv",
  NEEDS_EEJ: "t-eej",
  NO_SINGLE_UNIQUE_JUNCTION: "t-hard",
};
