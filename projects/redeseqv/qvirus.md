# MinION Mk1B

## Google Cloud

Install `Google Cloud`:

    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | sudo tee -a /etc/apt/sources.list.d/google-cloud-sdk.list
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | sudo apt-key --keyring /usr/share/keyrings/cloud.google.gpg add -
    sudo apt-get update && sudo apt-get install google-cloud-sdk

Test `Google Cloud` installation:

    gcloud

Download the secret keyfile (only laise.moraes@fiocruz.br):

    https://fiocruzbr.sharepoint.com/:u:/s/PVM-IGM-BackupMiseq/Ee7UfqQaV-5GgK1YQs8RmIcBR09nzd4XLv2WtAOiMkm1aQ?email=laise.moraes%40fiocruz.br&e=UJxtGR

Run gcloud auth login

    gcloud auth activate-service-account rede-qvirus-ba-salvador@viraseq.iam.gserviceaccount.com \
        --key-file=secret-salvador-ba.json \
        --no-user-output-enabled \
        --verbosity error

List files and folders with ls

    gsutil ls gs://rede-qvirus-ba-salvador-550374

Upload the sequencing directory

    gsutil cp -r <seq_dir> gs://rede-qvirus-ba-salvador-550374

Send e-mail to `suporte.qvirus@mendelics.com.br` and inform about the upload

## Mendelics Rede QVirus

Save metadata sequencing as `csv` and check:
- separator: change comma to semicolon
- metodo_coleta: change `Nasopharyngeal Swab` to `swab-nf`
- ct_1: change dot to comma
- data_coleta: `dd/mm/yyyy` format
- do not exclude empty columns
- do not add positive and negative controls

Login to `https://rede.qvirus.com.br/auth/login` and submit the sequencing metadata