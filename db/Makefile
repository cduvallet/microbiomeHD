
all: zenodo_info.html

# version_changes.md, zenodo_text.md and dataset_info.yaml are user-provided

zenodo_info.html: zenodo_text.md info_table.md version_changes.md
	pandoc zenodo_text.md -o zenodo_text.html
	pandoc info_table.md -o info_table.html
	pandoc version_changes.md -o version_changes.html
	cat zenodo_text.html info_table.html version_changes.html > zenodo_description_final.html

info_table.md: ../src/util/list_zenodo_datasets.py dataset_info.yaml
	python $< ~/keys/zenodo.pem dataset_info.yaml $@
