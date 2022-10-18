# PVMSEQ-DATABASES

## Win+R
```sh
cmd
```

## cmd.exe [DOS]
```sh
xcopy "N:\Grupos\DiagCOVID19\Sistemas\Soroteca\ControledeAmostras_FioCruz_be.mdb" "G:\" /y
```

```sh
xcopy "G:\ControledeAmostras_FioCruz_be.mdb" "D:\OneDrive - FIOCRUZ\Sequenciamento\BANCO_DE_DADOS\SOROTECA" /y
```

```sh
"D:\OneDrive - FIOCRUZ\Sequenciamento\SCRIPTS\PVM-SEQ_SEROTECA_D_DRIVE.pgm7"
```

## Win+R
```sh
wt
```

## wt.exe [WSL2]
```sh
PVMSEQ-DATABASES
```

# PVMSEQ-EXTRACTION_DATE

## Win+R
```sh
wt
```

## wt.exe [WSL2]
```sh
nano $HOME/bin/PVMSEQ-EXTRACTION_DATE
```

```sh
PVMSEQ-EXTRACTION_DATE
```

# PVMSEQ-DATA

## Win+R
```sh
wt
```

## wt.exe [WSL2]
```sh
cd $HOME/PVM_SEQ/CORRIDAS/DOCUMENTOS/IGM_PVM_LIBRARYyyyymmdd
```

```sh
PVMSEQ-DATA ../../SAMPLE_SHEETS/*.csv
```

# VIGEAS-ILLUMINA

## Win+R
```sh
wt
```

## wt.exe [WSL2]
```sh
UPDATE
vigeas-illumina -u
```

```sh
LIBRARY=IGM_PVM_LIBRARYyyyymmdd
bs download project --no-metadata --summary --extension=fastq.gz -o $HOME/BaseSpace/"$LIBRARY" -n "$LIBRARY"
bs download run --no-metadata --summary -o $HOME/BaseSpace/"$LIBRARY"_SAV -n "$LIBRARY"
```

```sh
LIBRARY=IGM_PVM_LIBRARYyyyymmdd
vigeas-illumina -w 1 -t 12 -s $HOME/PVM_SEQ/CORRIDAS/SAMPLE_SHEETS/"$LIBRARY".csv -i $HOME/BaseSpace/"$LIBRARY"
```

# PVMSEQ-REPORT

## PVM-SEQ_REDCap_IGM_PVM_LIBRARYyyyymmdd.xls >> Text (Tab delimited) (*.txt)

## Win+R
```sh
wt
```

## wt.exe [WSL2]
```sh
cd $HOME/PVM_SEQ/CORRIDAS/DOCUMENTOS/IGM_PVM_LIBRARYyyyymmdd
```

```sh
PVMSEQ-REPORT PVM-SEQ_REDCap_IGM_PVM_LIBRARY*.txt /home/lpmor22/IGM_SARSCOV2/IGM_PVM_LIBRARY*_depth10X_ANALYSIS/IGM_PVM_LIBRARY*.consensus.*.fasta
```

```sh
Seguem métricas dos controles e mocks:
CP: 00.00%
CP2: 00.00%
CP3: 00.00%
CnCDNA: 00.00%
CnCDNA2: 00.00%
CnCDNA3: 00.00%
CnPCR: 00.00%
CnPCR2: 00.00%
CnPCR3: 00.00%
MOCK: 00.00%
MOCK2: 00.00%
MOCK3: 00.00%
```

```sh
Planilha atualizada com métricas de qualidade da corrida
```

```sh
000 sequências submetidas no GISAID
000 sequências submetidas no relatório da Rede Genômica
```
