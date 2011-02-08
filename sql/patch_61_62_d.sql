# patch_61_62_s.sql
#
# title:  add experimental_group.description
#
# description:
# Add description field to experimental_group table

ALTER table experimental_group ADD `description` varchar(255) default NULL;

# patch identifier
INSERT INTO meta (species_id, meta_key, meta_value) VALUES (NULL, 'patch', 'patch_61_62_d.sql|experimental_group.description');


