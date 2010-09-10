# create the hapmap project

mk.project('HapMap', 'International HapMap Project');

# load Phase 2 and 3 HapMap subject data

data(hapmap)
mk.subject('HapMap', hapmap.subjects)

mk.subject.attr('HapMap', hapmap.subjects[-1],
                c('Family Identifier',
                  'Subject Name of Father',
                  'Subject Name of Mother',
                  'Reference Gender',
                  'Plate in Phase 2 HapMap',
                  'Sample Population Identifier'))
store.subject.data('HapMap', hapmap.subjects)

