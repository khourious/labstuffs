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

    xcopy "N:\Grupos\DiagCOVID19\Sistemas\Soroteca\ControledeAmostras_FioCruz_be.mdb" "G:\" /y

.. code:: bat

    xcopy "G:\ControledeAmostras_FioCruz_be.mdb" "C:\OneDrive\OneDrive - FIOCRUZ\Sequenciamento\BANCO_DE_DADOS\SOROTECA" /y

.. code:: bat

    xcopy "G:\ControledeAmostras_FioCruz_be.mdb" "D:\OneDrive - FIOCRUZ\Sequenciamento\BANCO_DE_DADOS\SOROTECA" /y

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

IGM_SARSCOV2
============

=====
Win+R
=====
.. code:: bat

    wt

=============
wt.exe [WSL2]
=============
.. code:: bash

    bs download project --no-metadata --summary --extension=fastq.gz -o $HOME/BaseSpace/IGM_PVM_LIBRARYyyyymmdd -n IGM_PVM_LIBRARYyyyymmdd

.. code:: bash

    UPDATE

.. code:: bash

    igm_sarscov2 -u

.. code:: bash

    bs download run --no-metadata --summary -o $HOME/BaseSpace/IGM_PVM_LIBRARYyyyymmdd_SAV -n IGM_PVM_LIBRARYyyyymmdd

.. code:: bash

    igm_sarscov2 -w 1 -t 16 -p ARTIC_V4-1 -i $HOME/BaseSpace/IGM_PVM_LIBRARYyyyymmdd 

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

    PVMSEQ-REPORT PVM-SEQ_REDCap_IGM_PVM_LIBRARY*.txt /home/lpmor22/IGM_SARSCOV2/IGM_PVM_LIBRARY*_depth10X_ANALYSIS/IGM_PVM_LIBRARY*.consensus.*.fasta


.. code:: bash
.. code:: bash
.. code:: bash
.. code:: bash
.. code:: bash
.. code:: bash
.. code:: bash
.. code:: bash
.. code:: bash
.. code:: bash
.. code:: bash
.. code:: bash
