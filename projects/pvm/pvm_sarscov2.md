# Plataforma de Vigilância Molecular

## Análises do sequenciamento de genoma completo de SARS-CoV-2

### Requisitos do sistema

|     |                                      |
| --- | ------------------------------------ |
| CPU | 8ª geração Intel Core i5 ou superior |
| RAM | 8 GB ou superior                     |
| SO  | Windows 10 (64-bit) ou superior      |
| WSL | versão 2                             |

### Programas necessários

- [Epi Info](https://www.cdc.gov/epiinfo/support/por/pt_downloads.html)
- [Illumina Sequencing Analysis Viewer (SAV)](https://support.illumina.com/sequencing/sequencing_software/sequencing_analysis_viewer_sav/downloads.html)
- [Microsoft Office 365 Educação](https://www.microsoft.com/pt-br/education/products/office)
- [Microsoft OneDrive](https://www.microsoft.com/pt-br/microsoft-365/onedrive/download)
- [Microsoft Teams](https://www.microsoft.com/pt-br/microsoft-teams/download-app)
- [Notepadd++](https://notepad-plus-plus.org/downloads/)
- [VIral GEnome ASsembly pipelines for WGS (vigeas)](https://github.com/khourious/vigeas)
- [Windows Subsystem for Linux 2 (WSL2)](https://learn.microsoft.com/pt-br/windows/wsl/install)

### Configuração do Linux para as análises e montagem dos genomas e relatórios

Logar no MS Teams e criar atalho dentro do OneDrive para o diretório de sequenciamento da PVM. Após sincronizar o OneDrive, checar se o diretório `Sequenciamento` está presente.

```
MS Teams -> PVM-IGM -> Sequenciamento -> Arquivos -> Adicionar atalho ao OneDrive
```

No **WSL2**, realizar os seguintes procedimentos abaixo:

```sh
ln -s /mnt/c/OneDrive/OneDrive\ -\ FIOCRUZ/ OneDrive # criar atalho para o diretório do OneDrive
ln -s /mnt/c/OneDrive/OneDrive\ -\ FIOCRUZ/Sequenciamento PVM_SEQ # criar atalho para o diretório Sequenciamento presente no OneDrive
ln -s /mnt/c/OneDrive/OneDrive\ -\ FIOCRUZ/Sequenciamento\ Backup/ PVM_SEQ_BKP # criar atalho para o diretório Sequenciamento Backup presente no OneDrive
[ ! -d $HOME/bin ] && mkdir $HOME/bin # criar diretório bin no $HOME do usuário do Linux
ln -s $HOME/OneDrive/Sequenciamento/SCRIPTS/* $HOME/bin/ # criar atalho para os scripts dos relatórios dento do $HOME/bin
source $HOME/.$(ps -p $$ -ocomm=)rc # recarregar o perfil de configuração do shell
[ ! -d /mnt/c/BaseSpace/ ] && mkdir /mnt/c/BaseSpace/ # criar diretório para armazenar os dados baixados do BaseSpace
ln -s /mnt/c/BaseSpace/ BaseSpace # criar atalho de acesso do diretório BaseSpace no $HOME do usuário do Linux
git clone --recursive https://github.com/khourious/vigeas.git # clonar o repositório do pipeline de montagem dos genomas
cd vigeas # entrar no diretório vigeas
chmod 700 -R INSTALL # dar permissões completas ao arquivo de instalação das dependências do vigeas 
bash INSTALL # rodar a instalação do vigeas
```

### Atualização das bases de dados utilizadas para os relatórios

É necessario copiar o arquivo `ControledeAmostras_FioCruz_be.mdb` localizado na intranet da FIOCRUZ `\IGM-FS\Arquivos\Grupos\DiagCOVID19\Sistemas\Soroteca` para o diretório `\OneDrive\OneDrive - FIOCRUZ\Sequenciamento\BANCO_DE_DADOS\SOROTECA`.

Futuramente, o arquivo `ControledeAmostras_FioCruz_be.mdb` estará disponível dentro do MS Teams e estará localizado em `\OneDrive\OneDrive - FIOCRUZ\Soroteca`.

Após copiar o `ControledeAmostras_FioCruz_be.mdb`, abrir o `Executar` utiizando o atalho `Windows + R` e rodar o `Prompt de Comando do Windows`:

```
cmd
```

No `cmd`: rodar script do Epi Info para exportar os dados do Sistema de Soroteca da PVM:

```
"C:\OneDrive\OneDrive - FIOCRUZ\Sequenciamento\SCRIPTS\PVM-SEQ_SEROTECA_C_DRIVE.pgm7"
```

No `WSL2`: rodar script para exportar e organizar todas as bases de dados utilizadas para os relatórios:

```sh
PVMSEQ-DATABASES
```

### Requisição da lista de amostras para extração e sequenciamento do SARS-CoV-2

No `WSL2`: editar o campo de data em `$75>="yyyy-mm-dd"`, o qual yyyy-mm-dd representa a data mínima da busca das amostras. Salvar arquivo após edição utilizando `Ctrl + X -> Y -> Enter`:

```sh
nano $HOME/bin/PVMSEQ-EXTRACTION_DATE
```

No `WSL2`: rodar script para gerar a planilha de requisição de extração de amostras. O arquivo será salvo em `\OneDrive\Sequenciamento\REQUISICOES_SEQ`:

```sh
PVMSEQ-EXTRACTION_DATE
```

### Montagem do genomas de SARS-CoV2

No `WSL2`: identificar o nome da biblioteca em `LIBRARY=IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd` e abrir arquivo da samplesheet do sequenciamento:

```sh
LIBRARY=IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd # criar array com o nome da biblioteca de sequenciamento
dos2unix $HOME/PVM_SEQ/CORRIDAS/SAMPLE_SHEETS/"$LIBRARY".csv
nano $HOME/PVM_SEQ/CORRIDAS/SAMPLE_SHEETS/"$LIBRARY".csv # editar samplesheet da biblioteca de sequenciamento
```

Avaliar a samplesheet do sequenciamento de acordo com os seguintes critérios:
- As amostras e controles não podem conter `-` ou `_` uma vez que estes caracteres são utilizados pelo script de montagem como delimitadores de arquivos
- As amostras devem ser identificadas pelo tracking ID biobanco
- Identificador dos contoles devem sempre conter caractere numérico (*i.e.* MOCK1, CNCDNA1, CNPCR1, CP1)
- A coluna descrição deve conter a informação do esquema de primer utilizado (*i.e.* ARTIC_V4_1).

No `WSL2`:

```sh
sudo apt-get -y update # atualizar lista de pacotes do linux
sudo apt-get -y full-upgrade # atualizar o linux e dependências instaladas
sudo apt-get autoremove # remover dependências que não são mais necessárias
sudo apt-get auto-clean # remover arquivos de instalações de dependências
sudo apt-get -y purge $(dpkg -l | awk '/^rc/ {print $2}') # remover arquivos de instalações de dependências que o auto-clean não consegue resolver
sudo apt-get check # checar se há dependências quebradas
conda clean -ay # limpar o cachê do conda
vigeas-illumina -u # atualizar as dependências utililizadas pelos ambientes do vigeas-illumina
```

No `WSL2`:

```sh
bs download project --no-metadata --summary --extension=fastq.gz -o $HOME/BaseSpace/"$LIBRARY" -n "$LIBRARY" # baixar os arquivos fastQ
bs download run --no-metadata --summary -o $HOME/BaseSpace/"$LIBRARY"_SAV -n "$LIBRARY" # baixar os arquivos de qualidade da corrida
```

No `WSL2`:

```sh
vigeas-illumina -w 1 -s $HOME/PVM_SEQ/CORRIDAS/SAMPLE_SHEETS/"$LIBRARY".csv -i $HOME/BaseSpace/"$LIBRARY" -d 10 # rodar o vigeas para realizar a montagem dos genomas de SARS-CoV-2
```

### Relatórios da PVM -- gerar arquivos para montagem do relatório REDCap

No `WSL2`:

```sh
cd $HOME/PVM_SEQ/CORRIDAS/DOCUMENTOS/"$LIBRARY" # entrar no diretório dos documentos da biblioteca de sequenciamento
for f in *\ *; do mv "$f" "${f// /_}"; done # renomear os arquivos para não apresentarem espaços
PVMSEQ-DATA $HOME/PVM_SEQ/CORRIDAS/SAMPLE_SHEETS/"$LIBRARY".csv # rodar script para gerar os documentos que serão utilizados na montagem dos relatórios
```

Os arquivos gerados pelo script são:
- `IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd_ASSEMBLY`: contém esquema de sequenciamento e versão dos programas utilizados para montagem dos genomas
- `IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd_BIOBANCO`: contém o tracking ID biobanco das amostras
- `IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd_PRIMERS`: contém o tracking ID biobanco e esquema de primers
- `IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd_REC`: contém o número REC equivalente para cada tracking ID biobanco
- `PVM-SEQ_REDCap_IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd`: relatório REDCap
- `SolicitacaoMDB_ViaBiobanco_yyyy-mm-dd`: arquivo auxiliar para montagem do relatório REDCap
- `SolicitacaoGal29_ViaGal_yyyy-mm-dd`: arquivo auxiliar para montagem do relatório REDCap

### Relatórios da PVM -- editar relatório REDCap

Abrir a planilha `PVM-SEQ_REDCap_IGM_PVM_LIBRARYyyyymmdd.tsv` utilizando o `MS Excel` e manipular o arquivo de acordo com os seguintes critérios:
- Transformar a disposição dos dados das colunas `gal_id` e `cns` em número e retirar os números decimais adicionados
- Copiar os dados para `notepadd++`
- Transformar a disposição dos dados da planilha em `texto`
- Colar os dados do `notepadd++` de volta para a planilha
- Salvar a planilha `PVM-SEQ_REDCap_IGM_PVM_LIBRARYyyyymmdd.tsv` como `PVM-SEQ_REDCap_IGM_PVM_LIBRARYyyyymmdd.xls`
- Completar os dados faltantes:
 - **seq_id**: confirmar ID no REDCap utilizando o arquivo `REDCap_SequenciamentoDeSARS_DATA_yyyy-mm-dd` disponível no diretório `\OneDrive\Sequenciamento\BANCO_DE_DADOS` e seguir a númeração para as novas entradas
 - **gal_id**: utilizar as planilhas `SolicitacaoGal29_ViaGal_yyyy-mm-dd.csv` e `SolicitacaoMDB_ViaBiobanco_yyyy-mm-dd` para completar os dados faltantes
 - **req_sequenc**: adicionar a sigla do requisitante do sequenciamento
   - `HSR`: Hospital São Rafael
   - `LABCOV`: Laboratório de Diagnóstico Molecular da COVID-19 do CCS/UFRB
   - `LACEN-BA`: Laboratório Central da Bahia
   - `LJC`: Laboratório Jaime Cerqueira
   - `PVM`: Plataforma de Vigilância Molecular do IGM/FIOCRUZ
 - **biobanco_seq / primer_schem**: utilizar o arquivo `IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd_PRIMERS` para obter as informações de biobanco e esquema de primers
 - **pipe_assembly / depth_assembly**: utilizar o arquivo `IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd_ASSEMBLY` para obter as informações de sequenciamento, montagem e profundidade utilizada para geração do genoma consenso
 - **dt_coleta**: dispor a data de coleta no formato yyyy-mm-dd