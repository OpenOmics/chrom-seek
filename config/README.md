
## Notes for **genome.json** file

- `EFFECTIVEGENOMESIZE` for each reference (hg38, hg19, and mm10) can be found here [deepTools](https://deeptools.readthedocs.io/en/develop/content/feature/effectiveGenomeSize.html).

- `HOMER_REF` is a path to pre-compiled homer references for particular genomes

> [!WARNING]  
> If you are editing genomes.json, particularly `HOMER_REF` or `ALIAS`, please ensure the `ALIAS` alias you select exists
> at `/fdb/homer/genomes/$ALIAS` or update the directory to the homer preparsed files.