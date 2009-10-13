# create the demo project

mk.project('Demo','Demo Project')

# load descriptions of 270 Phase 2 HapMap subjects

data(hapmap)
demo.subjects <- subset(hapmap.subjects, !is.na(plate))
demo.subjects$panel <- factor(demo.subjects$panel)
mk.subject('Demo', demo.subjects)

# set up HapMap subject attributes

mk.subject.attr('Demo', demo.subjects[-1],
                c('Family Identifier',
                  'Subject Name of Father',
                  'Subject Name of Mother',
                  'Reference Gender',
                  'Plate in Phase 2 HapMap',
                  'Sample Population Identifier'))
store.subject.data('Demo', demo.subjects)

#-----------------------------------------------------------

# first demo dataset: 5K SNPs distributed across the genome

# describe the platform and dataset

mk.platform('Demo_Set_1', 'Demo 10K SNP Set')
mk.assay.group('Demo_Set_1', 'Des01', 'Array 1 of 4')
mk.assay.group('Demo_Set_1', 'Des02', 'Array 2 of 4')
mk.assay.group('Demo_Set_1', 'Des03', 'Array 3 of 4')
mk.assay.group('Demo_Set_1', 'Des04', 'Array 4 of 4')
mk.mapping('Demo_Set_1', 'Demo_Set_1_b36',
           'Platform 1 Mapping to NCBI Build 36', 'ncbi_b36')
mk.dataset('Demo_1', 'Demo', 'Demo_Set_1', 'Demo (HapMap) Dataset #1')
mk.flag('assay_data', 'Demo_1', 1, 'pass', 'Passed QC filters')

# now load the data

data(demo_01)
mk.sample('Demo_1', samples.01)
mk.assay('Demo_Set_1', assay.def.01)
mk.assay.position('Demo_Set_1', , assay.pos.01)
mk.assay.data('Demo_1', assay.dat.01)

#-----------------------------------------------------------

# second dataset: ~400 SNPs from small interval on chr21

# describe the platform and dataset

mk.platform('Demo_Set_2', 'Demo Chr21 SNPs')
mk.assay.group('Demo_Set_2', 'Des11', 'Array 11')
mk.mapping('Demo_Set_2', 'Demo_Set_2_b36',
           'Platform 2 Mapping to NCBI Build 36', 'ncbi_b36')
mk.dataset('Demo_2', 'Demo', 'Demo_Set_2', 'Demo (HapMap) Dataset #2')
mk.flag('assay_data', 'Demo_2', 1, 'pass', 'Passed QC filters')

# now load the data: note that there are some duplicate samples

data(demo_02)
mk.sample('Demo_2', samples.02)
mk.sample.attr('Demo_2', samples.02['is_dup'])
store.sample.data('Demo_2', samples.02[c('sample.name','is_dup')])

mk.assay('Demo_Set_2', assay.def.02)
mk.assay.position('Demo_Set_2', , assay.pos.02)
mk.assay.data('Demo_2', assay.dat.02)

