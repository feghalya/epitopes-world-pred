#!/usr/bin/env python3

import sys

from camaptools.MLP import Model as CAMAPModel


class Model(CAMAPModel):
    """camaptools Model with peptide discovery from FASTA files.

    discover_peptides() is taken from Codon_publication's model_tools.py (d1d7e34).
    """

    def discover_peptides(self, fasta_file, peptide_lengths, coding=True, whitelist=None):
        import pyGeno.tools.UsefulFunctions as uf
        import pyGeno.tools.parsers.FastaTools as ft

        fa = ft.FastaFile(fasta_file)

        for header, sequence_fa in fa:
            header = header.replace('>', '')
            # if there is a whitelise use it
            if whitelist is None or sum([x in header for x in whitelist]):
                if 'NP_057857' in header:  # (dangerous: NP_057857XX)
                    sys.stderr.write('Skipping HIV Nef protein (has a UAG that behaves as UGG coding for W)\n')
                    continue
                # calculate score
                for peplen in peptide_lengths:
                    assert len(sequence_fa) > peplen
                    # detect frame
                    frame = None
                    reverse = None
                    start_found = None
                    end_found = None
                    if coding:
                        for cur_frame, fseq in enumerate(uf.translateDNA_6Frames(sequence_fa)):
                            if '*' not in fseq[:-1]:
                                # save previous values
                                prev_frame = frame
                                prev_reverse = reverse
                                prev_start_found = start_found
                                prev_end_found = end_found
                                # get new values
                                frame = cur_frame % 3
                                reverse = cur_frame >= 3
                                start_found = fseq[0] == 'M'
                                end_found = fseq[-1] == '*'
                                # call ambiguous frames
                                if prev_frame is not None:
                                    if (prev_start_found and not start_found) or \
                                            (prev_end_found and not end_found):
                                        frame = prev_frame
                                        reverse = prev_reverse
                                        start_found = prev_start_found
                                        end_found = prev_end_found
                                    elif (not prev_start_found and start_found) or \
                                            (not prev_end_found and end_found):
                                        pass
                                    else:
                                        sys.stderr.write(
                                                'Frames: %d (reverse:%s), %d (reverse:%s) ' % (
                                                    prev_frame, str(prev_reverse),
                                                    frame, str(reverse)) + header + '\n'
                                                )
                        assert frame is not None
                        # setup sequences
                        if reverse:
                            if frame:
                                sequence_fa = uf.reverseComplement(sequence_fa[:-frame])
                            else:
                                sequence_fa = uf.reverseComplement(sequence_fa)
                        else:
                            sequence_fa = sequence_fa[frame:]
                        # remove stop codon
                        if uf.translateDNA(sequence_fa[-3:]) == '*':
                            sequence_fa = sequence_fa[:-3]
                        sequences = [sequence_fa]
                        end_origins = [1]
                    else:
                        sequences = [sequence_fa, uf.reverseComplement(sequence_fa)]
                        end_origins = [1, -1]

                    for (sequence, end_origin) in zip(sequences, end_origins):
                        ntpeplen = peplen*3
                        for i in range(0, len(sequence) - ntpeplen, 3 if coding else 1):
                            nt_peptide = sequence[i:i+ntpeplen]
                            before = sequence[max(0, i-self.context):i]
                            after = sequence[min(len(sequence), i+ntpeplen):i+ntpeplen+self.context]
                            seq = before + nt_peptide + after

                            protein = uf.translateDNA(seq)
                            if '*' in protein:
                                if coding:
                                    sys.stderr.write(
                                            'Warning: Stop codon encountered in coding sequence %s\n' % header
                                            )
                            elif set(seq) - set('ATCGUatcgu'):
                                pass
                                sys.stderr.write('Non ATCGU encoutered in coding sequence\n')
                            else:
                                if len(before) < self.context:
                                    before = ''.join(['0'] * (self.context - len(before))) + before
                                assert len(before) == self.context

                                if len(after) < self.context:
                                    after = after + ''.join(['0'] * (self.context - len(after)))
                                assert len(after) == self.context

                                seq = before + nt_peptide + after

                                peptide = uf.translateDNA(nt_peptide)
                                score = float(self.get_seq_score_from_nt(seq)[0][1])
                                one_base_peptide_start = (i + 1) * end_origin

                                yield header, one_base_peptide_start, peptide, score
                            # break  # good spot for debugging
