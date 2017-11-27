import re

from smart_open import smart_open

__all__ = [
    'Interval',
    'RefGene',
]


class Interval(object):
    def __init__(
        self,
        chrom,
        start,
        end,
        strand='.',
        name='.',
        score=500,
        **metadata
    ):
        assert isinstance(start, int), 'Loci must be integers'
        assert isinstance(end, int), 'Loci must be integers'
        assert isinstance(name, str), 'Name must be a string'
        assert isinstance(score, (int, float)), 'Score must be a number'
        assert end - start > 0, 'Exclusive end must be greater than start'
        assert strand in ('.', '+', '-'), 'Strand must be ".", "+", "-" only'

        self.chrom = str(chrom)
        self.start = start
        self.end = end
        self.strand = strand
        self.name = name
        self.score = score
        self.metadata = metadata

    @property
    def SAM_interval(self):
        return f'{self.chrom}:{self.start}-{self.end}'

    @property
    def BED_interval(self):
        return '\t'.join(map(str, [
            self.chrom,
            self.start,
            self.end,
            self.name or '.',
            self.score,
            self.strand]))

    def get(self, item, default=None):
        temp = self.__dict__.copy()
        temp.update(self.metadata)
        return temp.get(item, default)

    def __getitem__(self, item):
        temp = self.__dict__.copy()
        temp.update(self.metadata)
        return temp.get(item, None)

    def __len__(self):
        return self.end - self.start

    def __eq__(self, other):
        return all((
            self.chrom == other.chrom,
            self.start == other.start,
            self.end == other.end,
            self.strand == other.strand))

    def __lt__(self, other):
        return all((
            self.chrom == other.chrom,
            self.start < other.start))

    def __le__(self, other):
        return all((
            self.chrom == other.chrom,
            self.start <= other.start))

    def __gt__(self, other):
        return all((
            self.chrom == other.chrom,
            self.end > other.end))

    def __ge__(self, other):
        return all((
            self.chrom == other.chrom,
            self.end >= other.end))

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'"{self.chrom}", '
            f'{self.start}, '
            f'{self.end}, '
            f'"{self.strand}")')

    def __str__(self):
        return self.__repr__()


class AminoRegion(Interval):
    def __init__(
        self,
        chrom,
        start,
        end,
        strand='.',
        name='.',
        score=500,
        **metadata
    ):
        super().__init__(
            chrom,
            start,
            end,
            strand=strand,
            name=name,
            metadata=metadata)


class CodingRegion(Interval):
    def __init__(
        self,
        chrom,
        start,
        end,
        strand='.',
        name='.',
        score=500,
        **metadata
    ):
        super().__init__(
            chrom,
            start,
            end,
            strand=strand,
            name=name,
            metadata=metadata)


class Exon(Interval):
    def __init__(
        self,
        chrom,
        start,
        end,
        strand='.',
        rank=None,
        frame_offset=-1,
        name='.',
        **metadata
    ):
        super().__init__(
            chrom,
            start,
            end,
            strand=strand,
            name=name,
            metadata=metadata)

        if (
            frame_offset is not None and
            not (isinstance(frame_offset, int) and
                 frame_offset in range(-1, 3))
        ):
            raise ValueError(
                'Frame offset must be None or integer and in set [-1, 2]')

        if not isinstance(rank, int) and rank > 0:
            raise ValueError('Rank must be a positive integer!')

        self.rank = rank
        self.frame_offset = None if frame_offset == -1 else frame_offset

    @property
    def is_UTR(self):
        return self.frame_offset is None

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'"{self.chrom}", '
            f'{self.start}, '
            f'{self.end}, '
            f'"{self.strand}", '
            f'rank={self.rank}, '
            f'frame_offset={self.frame_offset})')


class Transcript(Interval):
    def __init__(
        self,
        chrom,
        start,
        end,
        strand='.',
        name='.',
        accession=None,
        coding_start=None,
        coding_end=None,
        score=None,
        coding_start_status=None,
        coding_end_status=None,
        **metadata
    ):
        super().__init__(
            chrom,
            start,
            end,
            strand=strand,
            name=name,
            metadata=metadata)

        self.accession = accession
        self.transcript_start = start
        self.transcript_stop = end
        self.coding_start = coding_start
        self.coding_end = coding_end
        self.coding_start_status = coding_start_status
        self.coding_end_status = coding_end_status

        self._exons = []

    @property
    def num_exons(self):
        return len(self.exons)

    @property
    def exons(self):
        return self._exons

    @exons.getter
    def exons(self):
        return sorted(self._exons)

    @property
    def coding_intervals(self):
        intervals = []
        for exon in self.exons:
            if (
                exon.is_UTR or
                self.coding_start > exon.end or
                self.coding_end < exon.start
            ):
                continue

            if self.coding_start in range(exon.start, exon.end + 1):
                start = self.coding_start
                end = exon.end
            elif self.coding_end in range(exon.start, exon.end + 1):
                start = exon.start
                end = self.coding_end
            else:
                start = exon.start
                end = exon.end

            cds = CodingRegion(
                chrom=exon.chrom,
                start=start,
                end=end,
                strand=exon.strand)

            cds.exon = exon
            cds.transcript = self

            intervals.append(cds)

        return intervals

    @property
    def coding_length(self):
        return sum(len(cds) for cds in self.coding_intervals)

    def __repr__(self):
        return (
            f'{self.__class__.__name__}('
            f'"{self.chrom}", '
            f'{self.start}, '
            f'{self.end}, '
            f'"{self.strand}", '
            f'name="{self.name}", '
            f'accession="{self.accession}")')


class RefGene(object):
    def __init__(self, path):
        self.path = path
        self._handle = None

    @staticmethod
    def _line_to_transcript(line):
        (_, accession, chrom, strand, transcript_start, transcript_stop,
         coding_start, coding_end, num_exons, exon_starts, exon_ends,
         score, alt_name, coding_start_status, coding_end_status,
         exon_frames) = line

        transcript = Transcript(
            chrom,
            int(transcript_start),
            int(transcript_stop),
            strand,
            name=alt_name,
            accession=accession,
            coding_start=int(coding_start),
            coding_end=int(coding_end),
            score=score,
            coding_start_status=coding_start_status,
            coding_end_status=coding_end_status)

        exon_starts = exon_starts.split(',')
        exon_ends = exon_ends.split(',')
        exon_frames = exon_frames.split(',')

        if strand == '-':
            exon_ranks = range(int(num_exons), 0, -1)
        else:
            exon_ranks = range(1, int(num_exons) + 1)

        for start, end, frame_offset, rank in zip(
            exon_starts, exon_ends, exon_frames, exon_ranks
        ):
            if any(_ == '' for _ in (start, end, frame_offset)):
                continue

            exon = Exon(
                chrom=chrom,
                start=int(start),
                end=int(end),
                strand=strand,
                rank=rank,
                frame_offset=int(frame_offset))

            exon.transcript = transcript

            transcript._exons.append(exon)

        return transcript

    def transcript_by_accession(self, accession, flags=re.IGNORECASE):
        pattern = re.compile(accession, flags)

        for transcript in self:
            if pattern.match(transcript.accession):
                yield transcript

    def transcript_by_name(self, name, flags=re.IGNORECASE):
        pattern = re.compile(name, flags)

        for transcript in self:
            if pattern.match(transcript.name):
                yield transcript

    def __iter__(self):
        self._handle = smart_open(str(self.path))
        return self

    def __next__(self):
        line = next(self._handle).decode('utf-8').strip().split('\t')
        return self._line_to_transcript(line)

    def __repr__(self):
        return f'RefSeq("{self.path}")'
