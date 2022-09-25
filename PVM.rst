PVMSEQ-DATABASES
================

=====
Win+R
=====
.. code:: bat

    cmd

=============
cmd.exe [DOS]
=============
.. code:: bat

    xcopy "M:\DiagCOVID19\Sistemas\Soroteca\ControledeAmostras_FioCruz_be.mdb" "A:\" /y

.. code:: bat

    xcopy "A:\ControledeAmostras_FioCruz_be.mdb" "D:\OneDrive - FIOCRUZ\Sequenciamento\BANCO_DE_DADOS\SOROTECA" /y

.. code:: bat

    "C:\OneDrive\OneDrive - FIOCRUZ\Sequenciamento\SCRIPTS\PVM-SEQ_SEROTECA_C_DRIVE.pgm7"

.. code:: bat

    "D:\OneDrive - FIOCRUZ\Sequenciamento\SCRIPTS\PVM-SEQ_SEROTECA_D_DRIVE.pgm7"

=====
Win+R
=====
.. code:: bat

    wt

=============
wt.exe [WSL2]
=============
.. code:: bash

    PVMSEQ-DATABASES

PVMSEQ-EXTRACTION_DATE
======================

=====
Win+R
=====
.. code:: bat

    wt

=============
wt.exe [WSL2]
=============
.. code:: bash

    nano $HOME/bin/PVMSEQ-EXTRACTION_DATE

.. code:: bash

    PVMSEQ-EXTRACTION_DATE

PVMSEQ-DATA
===========

=====
Win+R
=====
.. code:: bat

    wt

=============
wt.exe [WSL2]
=============
.. code:: bash

    cd $HOME/PVM_SEQ/CORRIDAS/DOCUMENTOS/IGM_PVM_LIBRARYyyyymmdd

.. code:: bash

    PVMSEQ-DATA ../../SAMPLE_SHEETS/*.csv

VIGEAS-ILLUMINA
===============

=====
Win+R
=====
.. code:: bat

    wt

=============
wt.exe [WSL2]
=============
.. code:: bash

    UPDATE

.. code:: bash

    vigeas-illumina -u

.. code:: bash

    vigeas-illumina -w 1 -t 16 -s $HOME/PVM_SEQ/CORRIDAS/SAMPLE_SHEETS/"$LIBRARY".csv -i $HOME/BaseSpace/"$LIBRARY"
    
.. code:: bash

    LIBRARY=IGM_PVM_LIBRARYyyyymmdd
    bs download project --no-metadata --summary --extension=fastq.gz -o $HOME/BaseSpace/"$LIBRARY" -n "$LIBRARY"
    bs download run --no-metadata --summary -o $HOME/BaseSpace/"$LIBRARY"_SAV -n "$LIBRARY"

PVMSEQ-REPORT
=============

==========================================================================
PVM-SEQ_REDCap_IGM_PVM_LIBRARYyyyymmdd.xls >> Text (Tab delimited) (*.txt)
==========================================================================

=====
Win+R
=====
.. code:: bat

    wt

=============
wt.exe [WSL2]
=============
.. code:: bash

    cd $HOME/PVM_SEQ/CORRIDAS/DOCUMENTOS/IGM_PVM_LIBRARYyyyymmdd

.. code:: bash

    PVMSEQ-REPORT PVM-SEQ_REDCap_IGM_PVM_LIBRARY*.txt $HOME/vigeas/IGM_PVM_LIBRARY*_ANALYSIS/IGM_PVM_LIBRARY*.consensus.*.fasta

.. code:: bash

    Seguem métricas dos controles e mocks:
    # CP: 00.00%
    # CP2: 00.00%
    # CP3: 00.00%
    # CnCDNA: 00.00%
    # CnCDNA2: 00.00%
    # CnCDNA3: 00.00%
    # CnPCR: 00.00%
    # CnPCR2: 00.00%
    # CnPCR3: 00.00%
    # MOCK: 00.00%
    # MOCK2: 00.00%
    # MOCK3: 00.00%

.. code:: bash

    Planilha atualizada com métricas de qualidade da corrida

.. code:: bash

    000 sequências submetidas no GISAID
    000 sequências submetidas no relatório da Rede Genômica

.. code:: bash
.. code:: bash
.. code:: bash
.. code:: bash
.. code:: bash
.. code:: bash
.. code:: bash
.. code:: bash
.. code:: bash
