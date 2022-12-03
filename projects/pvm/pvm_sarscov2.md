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

A montagem dos genoma demora cerca de 5 minutos por genoma em um computador com *hardware* de 9ª geração Intel Core i7 com 16 GB de memória RAM.

Ao final da montagem dos genomas, avaliar as seguintes situações:
- `IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd.log.*. yyyy-mm-dd.txt`: se há erros na análise
- `IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd.consensus.*.yyyy-mm-dd.fasta`: se quantidade de genomas está em conformidade com a samplesheet e se estão com tamanho de 29.903 pb
- `IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd.coverage.*.yyyy-mm-dd.pdf`: se quantidade de plots de cobertura e profundidade está em conformidade com a samplesheet e se houve eventuais problemas de montagem
- `IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd.summary.SARS-CoV-2.*.yyyy-mm-dd.txt`: se quantidade de genomas está em conformidade com a samplesheet e se as métricas de montagem estão completas

No `WSL2`:

```sh
mv $HOME/vigeas/"$LIBRARY"_ANALYSIS $HOME/PVM_SEQ_BKP/ANALISES # mover a análise para o diretório de backup do sequenciamento
```

### Relatório REDCap

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

Abrir a planilha `PVM-SEQ_REDCap_IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd.tsv` utilizando o `MS Excel` e manipular o arquivo de acordo com os seguintes critérios:
- Transformar a disposição dos dados das colunas `gal_id` e `cns` em número e retirar os números decimais adicionados
- Copiar os dados para `notepadd++`
- Transformar a disposição dos dados da planilha em `texto`
- Colar os dados do `notepadd++` de volta para a planilha
- Salvar a planilha `PVM-SEQ_REDCap_IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd.tsv` em formato "Excel 97-2003 Workbook (*.xls)" com nome `PVM-SEQ_REDCap_IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd.xls`
- Completar os dados faltantes:
  - **seq_id**: confirmar ID no REDCap utilizando o arquivo `REDCap_SequenciamentoDeSARS_DATA_yyyy-mm-dd` disponível no diretório `\OneDrive\Sequenciamento\BANCO_DE_DADOS` e seguir a númeração para as novas entradas
  - **gal_id**: utilizar as planilhas `SolicitacaoGal29_ViaGal_yyyy-mm-dd.csv` e `SolicitacaoMDB_ViaBiobanco_yyyy-mm-dd` para completar os dados faltantes
  - **req_sequenc**: adicionar a sigla do requisitante do sequenciamento
    - `HSR`: Hospital São Rafael
    - `LABCOV`: Laboratório de Diagnóstico Molecular da COVID-19 do CCS/UFRB
    - `LACEN-BA`: Laboratório Central da Bahia
    - `LAPEM`: Laboratório de Patologia Estrutural e Molecular do IGM/FIOCRUZ
    - `LJC`: Laboratório Jaime Cerqueira
    - `PVM`: Plataforma de Vigilância Molecular do IGM/FIOCRUZ
  - **biobanco_seq / primer_schem**: utilizar o arquivo `IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd_PRIMERS` para obter as informações de biobanco e esquema de primers
  - **pipe_assembly / depth_assembly**: utilizar o arquivo `IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd_ASSEMBLY` para obter as informações de sequenciamento, montagem e profundidade utilizada para geração do genoma consenso
  - **dt_coleta**: dispor a data de coleta no formato yyyy-mm-dd
  - **gal_complete**: adicionar o valor 2 para todas as entradas
  - Abrir a planilha `IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd.summary.SARS-CoV-2.*.yyyy-mm-dd` utilizando o `MS Excel` e copiar as métricas de sequenciamento
  - Salvar as modificações e depois exportar a planilha em formato "Text (Tab delimited) (*.txt)" com nome `PVM-SEQ_REDCap_IGM_PVM_MISEQ_DNAP_LIBRARYyyyymmdd.txt`

No `WSL2`:

```sh
cd $HOME/PVM_SEQ/CORRIDAS/DOCUMENTOS/"$LIBRARY" # entrar no diretório dos documentos da biblioteca de sequenciamento
PVMSEQ-REPORT PVM-SEQ_REDCap_"$LIBRARY".txt $HOME/PVM_SEQ_BKP/ANALISES/"$LIBRARY"_ANALYSIS/"$LIBRARY".consensus.*.fasta # utilizar o relatório RECap e os arquivos fasta dos genomas para poder gerar os demais relatórios
```

Serão gerados os seguintes arquivos:
- `RelatorioCIEVS_yyyy-mm-dd.csv`: arquivo que será encaminhada como relatório para o **Centro de Informações Estratégicas em Vigilância em Saúde (CIEVS)** da Secretaria de Estado de Saúde da Bahia (SESAB). O arquivo está salvo no diretório `\OneDrive\Sequenciamento\RELATORIOS\CIEVS`.
- `hCoV-19_FIOCRUZ_BA_PVM_yyyymmdd.tsv` / `hCoV-19_FIOCRUZ_BA_PVM_yyyymmdd.fasta`: arquivos que serão utilizados para organizar a submissão dos genomas para o **GISAID**. Os arquivos estão salvos no diretório `\OneDrive\Sequenciamento\RELATORIOS\GISAID`.
- `RelatorioRedeGenomica_FIOCRUZBahia_yyyy-mm-dd.csv`: arquivo que será utilizado para organizar o relatório relatório da **Rede Genômica Fiocruz** da Fundação Oswaldo Cruz (FIOCRUZ). O arquivo está salvo no diretório `\OneDrive\Sequenciamento\RELATORIOS\REDE_GENOMICA`

### Submissão GISAID

Para a submissão são necessários um arquivo multifasta com os genomas e uma planilha com os metadados requeridos pelo GISAID. O modelo da planilha está disponível no diretório `\OneDrive\Sequenciamento\RELATORIOS`.
- Abrir o arquivo `hCoV-19_FIOCRUZ_BA_PVM_yyyymmdd.xlt`  utilizando o `MS Excel` e prencher com os dados presentes no arquivo `hCoV-19_FIOCRUZ_BA_PVM_yyyymmdd.tsv`
- Salvar a planilha `hCoV-19_FIOCRUZ_BA_PVM_yyyymmdd.xlt` em formato "Excel 97-2003 Workbook (*.xls)" com nome `hCoV-19_FIOCRUZ_BA_PVM_yyyymmdd.xls`
- Deletar `hCoV-19_FIOCRUZ_BA_PVM_yyyymmdd.tsv`

#### Submissão via web
- Entrar no endereço [https://www.epicov.org/epi3/frontend](https://www.epicov.org/epi3/frontend)
- Escrever login e senha -- **utilizar o login RKhour0**
- Acessar ambiente de submissão dos arquivos:
```
EpiCov -> Upload -> Batch Upload
```
- Adicionar em **Metadata as Excel or CSV** o arquivo `hCoV-19_FIOCRUZ_BA_PVM_yyyymmdd.xls`
- Adicionar em **Sequences as FASTA** o arquivo `hCoV-19_FIOCRUZ_BA_PVM_yyyymmdd.fasta`
- Modificar os arquivos em caso de alguma sequencia já ter sido adicionada previamente no GISAID

#### Submissão via cli3
- No `WSL2`: checar a validade do login via cli3.
```sh
cat $HOME/PVM_SEQ/RELATORIOS/GISAID/gisaid.authtoken | awk -F", " '{print $NF}' # abrir arquivo token de login e checar a data de expiração
```
- No `WLS2`: em caso de expiração de login, realizar a autenticação utilizando username RKhour0, senha e client-ID: cid-b3382c70dcc41:
```sh
cli3 authenticate
```
- No `WLS2`: converter o arquivo `hCoV-19_FIOCRUZ_BA_PVM_yyyymmdd.xls` para o formato *.csv:
```sh
ssconvert hCoV-19_FIOCRUZ_BA_PVM_$(date +'%Y%m%d').xls hCoV-19_FIOCRUZ_BA_PVM_$(date +'%Y%m%d').csv # converter arquivo *.xls para *.csv
```
- No `WLS2`: enviar os dados de submissão:
```sh
cli3 upload --database EpiCoV --metadata $HOME/PVM_SEQ/RELATORIOS/GISAID/hCoV-19_FIOCRUZ_BA_PVM_$(date +'%Y%m%d').csv --fasta $HOME/PVM_SEQ/RELATORIOS/GISAID/hCoV-19_FIOCRUZ_BA_PVM_$(date +'%Y%m%d').fasta
```
- Modificar os arquivos em caso de alguma sequencia já ter sido adicionada previamente no GISAID

**Somente prosseguir para os demais relatórios (CIEVS e Rede Genômica Fiocruz) após submissão dos genomas no GISAID**.

### Relatório CIEVS

Abrir o arquivo `RelatorioCIEVS_yyyy-mm-dd.csv` utilizando o `MS Excel` e avaliar de acordo com os seguintes critérios:
- Checar se a coluna `NUMERO DE REQUISICAO` está devidamente preenchida
- Checar se a coluna `VOC` contém a informação escrita da VOC em questão


#### E-mail para o CIEVS

Enviar e-mail para o CIEVS com o arquivo `RelatorioCIEVS_yyyy-mm-dd.csv` em anexo, com as seguintes informações:
- Destinatários
```
cievs.notifica@saude.ba.gov.br; pvm@fiocruz.br; ricardo_khouri@hotmail.com
```
- Assunto
```
Relatório CIEVS de linhagens de SARS-CoV-2 no estado da Bahia
```
- Corpo do e-mail
```
 Prezados,
 
 Segue anexo o relatório com a linhagem dos genomas de SARS-CoV-2 recuperado de amostras no estado da Bahia. Estamos a disposição para maiores esclarecimentos. Muito obrigado pela colaboração na rede de vigilância genômica no país.
```

### Relatório Rede Genômica Fiocruz

Para a submissão são necessários um arquivo *.pdf do relatório e uma planilha com os dados que compõem o relatório. O modelo do relatório está disponível no diretório `\OneDrive\Sequenciamento\RELATORIOS`.
- Abrir o arquivo modelo do relatório utilizando o `MS PowerPoint`. Utilizar os arquivos de acordo com o requisitante do sequenciamento:
  - **HSR**: `RelatorioRedeGenomica_FIOCRUZBahia_HSR.potx`
  - **LABCOV**: `RelatorioRedeGenomica_FIOCRUZBahia_LABCOV.potx`
  - **LACEN-BA** / **PVM**: `RelatorioRedeGenomica_FIOCRUZBahia_GAL.potx`
  - **LAPEM**: `RelatorioRedeGenomica_FIOCRUZBahia_LAPEM.potx`
  - **LJC**: `RelatorioRedeGenomica_FIOCRUZBahia_LJC.potx`
  - **PVM**: `RelatorioRedeGenomica_FIOCRUZBahia_GAL.potx`
- Adicionar as seguintes informações adicionais:
  - A data do dia da confecção do relatório é informada como `Versão`
  - Especificar o laboratório de origem no caso do `LACEN-BA` ou `PVM-IGM`
  - Substituir `XX` pela quantidade de genomas incluídos no relatório em questão
  - Substituir a versao do pangolin e da base de dados pela utilizada na montagem (*i.e.* Pangolin v.4.0.6 / PUSHER-v1.9)
  - Abrir o arquivo `RelatorioRedeGenomica_FIOCRUZBahia_yyyy-mm-dd.csv` e configurar os dados em: **Arial, fonte tamanho 10, centralizado**
  - Prencher as tabelas com os dados configurados da planilha `RelatorioRedeGenomica_FIOCRUZBahia_yyyy-mm-dd.csv`
  - Salvar o arquivo em formato `PowerPoint Presentation (*.pptx)` com nome `RelatorioRedeGenomica_FIOCRUZBahia_yyyy-mm-dd.pptx`
  - Exportar o arquivo em formato `PDF (*.pdf)` com nome `RelatorioRedeGenomica_FIOCRUZBahia_yyyy-mm-dd.pdf`
  - Na planilha `Produção RGF enviados para o ministério - fechamento toda a quinta-feira`, na aba `IGM`, adicionar a quantidade de genomas submetidos no **GISAID** e nos relatórios da Rede Genômica Fiocruz. A planilha está disponível no [Google Sheets](https://docs.google.com/spreadsheets/d/1UXb4GcDQ7iKnM92Z20iJaxWWgyw536H0By7RdrsqkRc/edit?pli=1#gid=0)

#### E-mail para a Rede Genômica Fiocruz # LACEN-BA / PVM-IGM

Enviar e-mail para os integrantes da Rede Genômica Fiocruz com os arquivos `RelatorioRedeGenomica_FIOCRUZBahia_yyyy-mm-dd.csv` e `RelatorioRedeGenomica_FIOCRUZBahia_yyyy-mm-dd.pdf` em anexo, com as seguintes informações:

- Destinatários
```
arnaldo.medeiros@saude.gov.br; camila.indiani@fiocruz.br; cievs.notifica@saude.ba.gov.br; elisa.cavalcante@ioc.fiocruz.br; fcm@ioc.fiocruz.br; greice.madeleine@saude.gov.br; gripe@saude.gov.br; lacen.clavep@saude.ba.gov.br; marilda.goncalves@fiocruz.br; miriam.livorati@saude.gov.br; mmsiq@ioc.fiocruz.br; notificasalvador@gmail.com; paola@ioc.fiocruz.br; pvm@fiocruz.br; ricardo.riccio@fiocruz.br; ricardo_khouri@hotmail.com; thiago.guedes@saude.gov.br; tiago.graf@fiocruz.br; vitor.martins@ioc.fiocruz.br; walquiria.almeida@saude.gov.br
```
- Assunto
```
Relatório de linhagens de SARS-CoV-2 no estado da Bahia
```
- Corpo do e-mail
```
 Prezados,
  
 Segue anexo o relatório com a linhagem dos genomas de SARS-CoV-2 recuperado de amostras no estado da Bahia. Estamos a disposição para maiores esclarecimentos. Muito obrigado pela colaboração na rede de vigilância genômica no país.
```

#### E-mail para a Rede Genômica Fiocruz # HSR

Enviar e-mail para os integrantes da Rede Genômica Fiocruz com os arquivos `RelatorioRedeGenomica_FIOCRUZBahia_yyyy-mm-dd_HSR.csv` e `RelatorioRedeGenomica_FIOCRUZBahia_yyyy-mm-dd_HSR.pdf` em anexo, com as seguintes informações:

- Destinatários
```
arnaldo.medeiros@saude.gov.br; brunosolanosouza@gmail.com; camila.bacia@hsr.com.br; camila.indiani@fiocruz.br; carolina.nonaka@hsr.com.br; cievs.notifica@saude.ba.gov.br; elisa.cavalcante@ioc.fiocruz.br; fcm@ioc.fiocruz.br; greice.madeleine@saude.gov.br; gripe@saude.gov.br; lacen.clavep@saude.ba.gov.br; marilda.goncalves@fiocruz.br; miriam.livorati@saude.gov.br; mmsiq@ioc.fiocruz.br; notificasalvador@gmail.com; paola@ioc.fiocruz.br; pvm@fiocruz.br; ricardo.riccio@fiocruz.br; ricardo_khouri@hotmail.com; thiago.guedes@saude.gov.br; tiago.graf@fiocruz.br; vitor.martins@ioc.fiocruz.br; walquiria.almeida@saude.gov.br
```
- Assunto
```
Relatório de linhagens de SARS-CoV-2 no estado da Bahia
```
- Corpo do e-mail
```
 Prezados,
  
 Segue anexo o relatório com a linhagem dos genomas de SARS-CoV-2 recuperados de amostras Hospital São Rafael. Para maiores esclarecimentos sobre as amostras, por favor, contactar Dra Camila Bacia (camila.bacia@hsr.com.br), infectologista do Hospital São Rafael. Muito obrigado pela colaboração na rede de vigilância genômica no país.
```