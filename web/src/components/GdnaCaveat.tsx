import { Alert } from "./icons";

/**
 * What a same-exon product costs, said where the pair is offered.
 *
 * Every other product on this site crosses an exon–exon junction, so contaminating genomic
 * DNA either fails to amplify or gives a visibly longer band. A product inside ONE exon has
 * no such guard: the same two sites sit uninterrupted in the genome and give the identical
 * band. It is offered only where no junction-crossing product can do the job — a
 * single-exon transcript, or a whole-transcript pair that has to include one — and never
 * without this beside it. Visible, not a tooltip: it changes what the user does at the bench.
 */
export default function GdnaCaveat({ children }: { children?: React.ReactNode }) {
  return (
    <p className="gdna-note" role="note">
      <span className="gdna-ic"><Alert /></span>
      <span>
        {children}{children ? " " : ""}
        This product lies <b>inside one exon</b>, so genomic DNA gives the same band.{" "}
        <b>DNase-treat the RNA</b> and run a <b>no-RT control</b> alongside.
      </span>
    </p>
  );
}
