--
-- Copyright (C) 2009, Perlegen Sciences, Inc.
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

--
-- Assay Definitions and Mapping Information
--

create table platform
(
  platform_id integer primary key,
  name varchar(64) unique not null,
  description varchar(255),
  created_by varchar(64),
  created_dt date
);

create table assay_group
(
  assay_group_id integer primary key,
  platform_id integer not null,
  name varchar(64) unique not null,
  description varchar(255),
  created_by varchar(64),
  created_dt date,
  foreign key (platform_id)
    references platform(platform_id) on delete cascade
);

-- create table assay_flag
-- (
--  platform_id integer not null,
--  position integer not null,
--  name varchar(64) not null,
--  description varchar(255),
--  created_by varchar(64),
--  created_dt date,
--  constraint af_pk primary key (platform_id, position),
--  constraint af_uniq unique (platform_id, name),
--  foreign key (platform_id)
--    references platform(platform_id) on delete cascade
-- );

create table assay
(
  assay_id integer primary key,
  platform_id integer not null,
  assay_group_id integer not null,
  name varchar(255) not null,
  -- flags integer,
  alleles varchar(255),
  probe_seq varchar(255),
  constraint assay_uniq unique (platform_id, name),
  foreign key (platform_id)
    references platform(platform_id) on delete cascade,
  foreign key (assay_group_id)
    references assay_group(assay_group_id) on delete cascade
);

create table mapping
(
  mapping_id integer primary key,
  platform_id integer not null,
  name varchar(64) unique not null,
  description varchar(255),
  assembly varchar(64),
  is_hidden integer,
  created_by varchar(64),
  created_dt date,
  foreign key (platform_id) references platform(platform_id)
);

-- create table assay_position_flag
-- (
--  mapping_id integer not null,
--  position integer not null,
--  name varchar(64) not null,
--  description varchar(255),
--  created_by varchar(64),
--  created_dt date,
--  constraint adf_pk primary key (mapping_id, position),
--  constraint adf_uniq unique (mapping_id, name),
--  foreign key (mapping_id)
--    references mapping(mapping_id) on delete cascade
-- );

create table assay_position
(
  mapping_id integer not null,
  assay_id integer not null,
  -- flags integer,
  scaffold varchar(64),
  position integer,
  strand char(1),
  ploidy char(1),
  dbsnp_rsid integer,
  dbsnp_orient char(1),
  constraint assay_pos_uniq unique (mapping_id, assay_id),
  foreign key (mapping_id)
    references mapping(mapping_id) on delete cascade,
  foreign key (assay_id) references assay(assay_id)
);

create index assay_position_idx on
assay_position(mapping_id, scaffold, position);

create index assay_position_rsid_idx on
assay_position(dbsnp_rsid);

--
-- Projects and Datasets
--

create table project
(
  project_id integer primary key,
  name varchar(64) unique not null,
  description varchar(255),
  is_hidden integer,
  created_by varchar(64),
  created_dt date
);

create table dataset
(
  dataset_id integer primary key,
  project_id integer not null,
  platform_id integer not null,
  name varchar(64) unique not null,
  description varchar(255),
  is_hidden integer,
  created_by varchar(64),
  created_dt date,
  foreign key (project_id) references project(project_id),
  foreign key (platform_id) references platform(platform_id)
);

--
-- Subjects, Samples and Attributes
--

create table subject
(
  subject_id integer primary key,
  project_id integer not null,
  name varchar(255) not null,
  constraint subject_uniq unique (project_id, name),
  foreign key (project_id) references project(project_id)
);

create table subject_attr
(
  subject_attr_id integer primary key,
  project_id integer not null,
  name varchar(64) not null,
  datatype varchar(64),
  levels varchar(1000),
  description varchar(255),
  is_hidden integer,
  created_by varchar(64),
  created_dt date,
  constraint subject_attr_uniq unique (project_id, name),
  foreign key (project_id) references project(project_id)
);

create table subject_value
(
  subject_id integer not null,
  subject_attr_id integer not null,
  value varchar(255),
  constraint subject_value_pk primary key (subject_id, subject_attr_id),
  foreign key (subject_id)
    references subject(subject_id) on delete cascade,
  foreign key (subject_attr_id)
    references subject_attr(subject_attr_id) on delete cascade
);

create table sample
(
  sample_id integer primary key,
  dataset_id integer not null,
  subject_id integer not null,
  name varchar(255) not null,
  gender char(1),
  position integer,
  constraint sample_uniq_1 unique (dataset_id, name),
  constraint sample_uniq_2 unique (dataset_id, position),
  foreign key (dataset_id) references dataset(dataset_id),
  foreign key (subject_id) references subject(subject_id)
);

create table sample_attr
(
  sample_attr_id integer primary key,
  dataset_id integer not null,
  name varchar(64) not null,
  datatype varchar(64),
  levels varchar(1000),
  description varchar(255),
  is_hidden integer,
  created_by varchar(64),
  created_dt date,
  constraint sample_attr_uniq unique (dataset_id, name),
  foreign key (dataset_id) references dataset(dataset_id)
);

create table sample_value
(
  sample_id integer not null,
  sample_attr_id integer not null,
  value varchar(255),
  constraint sample_value_pk primary key (sample_id, sample_attr_id),
  foreign key (sample_id)
    references sample(sample_id) on delete cascade,
  foreign key (sample_attr_id)
    references sample_attr(sample_attr_id) on delete cascade
);

--
-- Assay Genotype Data
--

create table assay_data_flag
(
  dataset_id integer not null,
  position integer not null,
  name varchar(64) not null,
  description varchar(255),
  created_by varchar(64),
  created_dt date,
  constraint adf_pk primary key (dataset_id, position),
  constraint adf_uniq unique (dataset_id, name),
  foreign key (dataset_id)
    references dataset(dataset_id) on delete cascade
);

create table assay_data
(
  assay_data_id integer primary key,
  dataset_id integer not null,
  assay_id integer not null,
  flags integer default 0,
  genotype clob,
  qscore blob,
  raw_data blob,
  constraint assay_data_uniq unique (dataset_id, assay_id),
  foreign key (dataset_id) references dataset(dataset_id),
  foreign key (assay_id) references assay(assay_id)
);

--
-- Association test results
--

create table test
(
  test_id integer primary key,
  dataset_id integer not null,
  name varchar(64) not null,
  description varchar(255),
  fit varchar(64),
  model varchar(255),
  term varchar(64),
  is_hidden integer,
  created_by varchar(64),
  created_dt date,
  constraint test_uniq unique (dataset_id, name, term),
  foreign key (dataset_id) references dataset(dataset_id)
);

create table test_result
(
  test_id integer,
  assay_data_id integer,
  pvalue float,
  effect float,
  stderr float,
  constraint test_result_pk primary key (test_id, assay_data_id),
  foreign key (test_id) references test(test_id) on delete cascade,
  foreign key (assay_data_id) references assay_data(assay_data_id)
);

--
-- Principal Components Analysis results
--

create table prcomp
(
  prcomp_id integer primary key,
  dataset_id integer not null,
  name varchar(64) not null,
  description varchar(255),
  fn_call varchar(255),
  components integer,
  samples integer,
  assays integer,
  is_hidden integer,
  created_by varchar(64),
  created_dt date,
  constraint prcomp_uniq unique (dataset_id, name),
  foreign key (dataset_id) references dataset(dataset_id)
);

create table prcomp_component
(
  prcomp_id integer not null,
  component integer not null,
  sdev float,
  foreign key (prcomp_id)
    references prcomp(prcomp_id) on delete cascade
);

create table prcomp_loading
(
  prcomp_id integer not null,
  component integer not null,
  sample_id integer not null,
  loading float,
  constraint psl_pk primary key (prcomp_id, component, sample_id),
  foreign key (prcomp_id)
    references prcomp(prcomp_id) on delete cascade,
  foreign key (sample_id) references sample(sample_id)
);
create sequence platform_id_seq;
create trigger insert_platform
before insert on platform for each row
begin
  select platform_id_seq.nextval
  into :new.platform_id from dual;
end;

create sequence assay_group_id_seq;
create trigger insert_assay_group
before insert on assay_group for each row
begin
  select assay_group_id_seq.nextval
  into :new.assay_group_id from dual;
end;

create sequence assay_id_seq;
create trigger insert_assay
before insert on assay for each row
begin
  select assay_id_seq.nextval
  into :new.assay_id from dual;
end;

create sequence mapping_id_seq;
create trigger insert_mapping
before insert on mapping for each row
begin
  select mapping_id_seq.nextval
  into :new.mapping_id from dual;
end;

create sequence project_id_seq;
create trigger insert_project
before insert on project for each row
begin
  select project_id_seq.nextval
  into :new.project_id from dual;
end;

create sequence dataset_id_seq;
create trigger insert_dataset
before insert on dataset for each row
begin
  select dataset_id_seq.nextval
  into :new.dataset_id from dual;
end;

create sequence subject_id_seq;
create trigger insert_subject
before insert on subject for each row
begin
  select subject_id_seq.nextval
  into :new.subject_id from dual;
end;

create sequence subject_attr_id_seq;
create trigger insert_subject_attr
before insert on subject_attr for each row
begin
  select subject_attr_id_seq.nextval
  into :new.subject_attr_id from dual;
end;

create sequence sample_id_seq;
create trigger insert_sample
before insert on sample for each row
begin
  select sample_id_seq.nextval
  into :new.sample_id from dual;
end;

create sequence sample_attr_id_seq;
create trigger insert_sample_attr
before insert on sample_attr for each row
begin
  select sample_attr_id_seq.nextval
  into :new.sample_attr_id from dual;
end;

create sequence assay_data_id_seq;
create trigger insert_assay_data
before insert on assay_data for each row
begin
  select assay_data_id_seq.nextval
  into :new.assay_data_id from dual;
end;

create sequence test_id_seq;
create trigger insert_test
before insert on test for each row
begin
  select test_id_seq.nextval
  into :new.test_id from dual;
end;

create sequence prcomp_id_seq;
create trigger insert_prcomp
before insert on prcomp for each row
begin
  select prcomp_id_seq.nextval
  into :new.prcomp_id from dual;
end;

