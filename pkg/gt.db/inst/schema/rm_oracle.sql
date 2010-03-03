--
-- Copyright (C) 2009, Perlegen Sciences, Inc.
-- Copyright (C) 2010, 23andMe, Inc.
-- 
-- Written by David A. Hinds <dhinds@sonic.net>
-- 
-- This is free software; you can redistribute it and/or modify it
-- under the terms of the GNU General Public License as published by
-- the Free Software Foundation; either version 3 of the license, or
-- (at your option) any later version.
-- 
-- This program is distributed in the hope that it will be useful,
-- but WITHOUT ANY WARRANTY; without even the implied warranty of
-- MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
-- GNU General Public License for more details.
-- 
-- You should have received a copy of the GNU General Public License
-- along with this program.  If not, see <http://www.gnu.org/licenses/>
-- 

drop trigger insert_prcomp;
drop sequence prcomp_id_seq;
drop trigger insert_test;
drop sequence test_id_seq;
drop trigger insert_assay_data;
drop sequence assay_data_id_seq;
drop trigger insert_sample_attr;
drop sequence sample_attr_id_seq;
drop trigger insert_sample;
drop sequence sample_id_seq;
drop trigger insert_subject_attr;
drop sequence subject_attr_id_seq;
drop trigger insert_subject;
drop sequence subject_id_seq;
drop trigger insert_dataset;
drop sequence dataset_id_seq;
drop trigger insert_project;
drop sequence project_id_seq;
drop trigger insert_mapping;
drop sequence mapping_id_seq;
drop trigger insert_assay;
drop sequence assay_id_seq;
drop trigger insert_platform;
drop sequence platform_id_seq;
drop table prcomp_loading;
drop table prcomp_component;
drop table prcomp;
drop table test_result;
drop table test;
drop table assay_data;
drop table assay_data_flag;
drop table sample_value;
drop table sample_attr;
drop table sample;
drop table subject_value;
drop table subject_attr;
drop table subject;
drop table dataset;
drop table project;
drop table assay_position;
drop table assay_position_flag;
drop table mapping;
drop table assay_flag;
drop table assay;
drop table platform;
drop table gtdb_option;
