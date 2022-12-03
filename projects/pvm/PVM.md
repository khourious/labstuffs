# Plataforma de Vigilância Molecular

## Análises do sequenciamento de genoma completo de SARS-CoV-2

### Requisitos do sistema

|     |                                      |
| --- | ------------------------------------ |
| CPU | 8ª geração Intel Core i5 ou superior |
| RAM | 8 GB ou superior                     |
| SO  | Windows 10 (64-bit) ou superior      |
| WSL | Windows Subsystem for Linux 2        |

### Programas necessários

- [Epi Info](https://www.cdc.gov/epiinfo/support/por/pt_downloads.html)
- [Illumina Sequencing Analysis Viewer (SAV)](https://support.illumina.com/sequencing/sequencing_software/sequencing_analysis_viewer_sav/downloads.html)
- [MS Teams](https://www.microsoft.com/pt-br/microsoft-teams/download-app)
- [MS OneDrive](https://www.microsoft.com/pt-br/microsoft-365/onedrive/download)
- [VIral GEnome ASsembly pipelines for WGS (vigeas)](https://github.com/khourious/vigeas)
- [WSL2](https://learn.microsoft.com/pt-br/windows/wsl/install)

### Configuração do Linux para as análises e montagem dos genomas e relatórios

Logar no MS Teams e criar atalho dentro do OneDrive para o diretório de sequenciamento da PVM. Após sincronizar o OneDrive, checar se o diretório `Sequenciamento` está presente.

```
MS Teams -> PVM-IGM -> Sequenciamento -> Arquivos -> Adicionar atalho ao OneDrive
```

No **WSL2**, realizar os seguintes procedimentos abaixo:

```sh
ln -s /mnt/c/OneDrive\ -\ FIOCRUZ/ OneDrive # criar atalho para o diretório do OneDrive
[ ! -d $HOME/bin ] && mkdir $HOME/bin # criar diretório bin no $HOME do usuário do Linux
ln -s $HOME/OneDrive/Sequenciamento/SCRIPTS/* $HOME/bin/ # criar atalho para os scripts dos relatórios dento do $HOME/bin
source $HOME/.$(ps -p $$ -ocomm=)rc # recarregar o perfil de configuração do shell
mkdir /mnt/c/BaseSpace/ # criar diretório para armazenar os dados baixados do BaseSpace
ln -s /mnt/c/BaseSpace/ BaseSpace # criar atalho de acesso do diretório BaseSpace no $HOME do usuário do Linux
git clone --recursive https://github.com/khourious/vigeas.git # clonar o repositório do pipeline de montagem dos genomas
cd vigeas # entrar no diretório vigeas
chmod 700 -R INSTALL # dar permissões completas ao arquivo de instalação das dependências do vigeas 
bash INSTALL # rodar a instalação do vigeas
```

### Atualização das bases de dados utilizadas para os relatórios

É necessario copiar o arquivo `ControledeAmostras_FioCruz_be.mdb` localizado na intranet da FIOCRUZ `\IGM-FS\Arquivos\Grupos\DiagCOVID19\Sistemas\Soroteca` para o diretório `\OneDrive\Sequenciamento\BANCO_DE_DADOS\SOROTECA`.

Futuramente, o arquivo `ControledeAmostras_FioCruz_be.mdb` estará disponível dentro do MS Teams e estará localizado em `\OneDrive\Diagnostico\Sistemas\Soroteca`.

Após copiar o `ControledeAmostras_FioCruz_be.mdb`, abrir o `Executar` utiizando o atalho `Windows + R` e rodar o `Prompt de Comando do Windows`:

```
cmd
```

No `cmd`: rodar script do Epi Info para exportar os dados do Sistema de Soroteca da PVM:

```
"C:\OneDrive - FIOCRUZ\Sequenciamento\SCRIPTS\PVM-SEQ_SEROTECA_D_DRIVE.pgm7"
```

No `WSL2`: rodar script para exportar e organizar todas as bases de dados utilizadas para os relatórios:

```sh
PVMSEQ-DATABASES
```

### Requisição da lista de amostras para extração e sequenciamento do SARS-CoV-2

No `WSL2`: editar o campo de data em `$75>="YYYY-MM-DD"`, o qual YYYY-MM-DD representa a data mínima da busca das amostras. Salvar arquivo após edição utilizando `Ctrl + X -> Y -> Enter`:

```sh
nano $HOME/bin/PVMSEQ-EXTRACTION_DATE
```

No `WSL2`: rodar script para gerar a planilha de requisição de extração de amostras. O arquivo será salvo em `\OneDrive\Sequenciamento\REQUISICOES_SEQ`:

```sh
PVMSEQ-EXTRACTION_DATE
```