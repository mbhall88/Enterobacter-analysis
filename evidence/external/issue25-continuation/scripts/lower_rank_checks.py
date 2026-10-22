"""Actual R2 omits genomovar fields; R1 may carry explicit null fields."""
def assert_no_genomovar(report):
 paths=[]
 for key in ['genomovar_id','diagnostic_would_assign_genomovar_id']:
  assert report.get(key) is None,'Populated genomovar field on ANIm rescue: '+key;paths.append({'path':key,'present':key in report,'value':report.get(key)})
 for section_name in ['genomovar','assignment']:
  section=report.get(section_name)
  if section is None:continue
  assert isinstance(section,dict),'Malformed optional lower-rank section'
  for key in ['genomovar_id','diagnostic_would_assign_genomovar_id']:
   assert section.get(key) is None,'Populated lower-rank field on ANIm rescue: '+section_name+'.'+key;paths.append({'path':section_name+'.'+key,'present':key in section,'value':section.get(key)})
 return paths
