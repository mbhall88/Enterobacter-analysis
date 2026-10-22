"""Pure budget/GPFS parsing and own-path filesystem accounting; no submissions."""
import os,pathlib
MIB=1024**2;GIB=1024**3

def parse_fileset_quota(text):
 line=next(x for x in text.splitlines() if x.split()[:2]==['scratch','FILESET']);left,right=line.split('|',1);b=left.split();f=right.split()
 return {'block_used_bytes':int(b[2])*1024,'block_soft_bytes':int(b[3])*1024,'block_hard_bytes':int(b[4])*1024,'block_in_doubt_bytes':int(b[5])*1024,'files_used':int(f[0]),'files_soft':int(f[1]),'files_hard':int(f[2]),'files_in_doubt':int(f[3])}

def own_usage(root):
 seen=set();size=0;files=0
 for parent,ds,fs in os.walk(root,followlinks=False):
  for name in ds+fs:
   p=pathlib.Path(parent)/name;s=p.lstat();files+=1
   if p.is_symlink():continue
   key=(s.st_dev,s.st_ino)
   if key not in seen:seen.add(key);size+=s.st_blocks*512
 return {'allocated_bytes':size,'file_and_directory_entries':files}

def admission(usage,quota,n,cap):
 assert 1<=n<=100
 duplicate_and_incoming=2*n*5*MIB;recovery=n*25*MIB;other=768*MIB;needed=duplicate_and_incoming+recovery+other
 inode_need=10000+30*n+2*n+1000
 return {'new_pairs':n,'own_usage':usage,'quota':quota,'cap_bytes':cap,'archive_and_collector_reserve_bytes':duplicate_and_incoming,'raw_recovery_reserve_bytes':recovery,'metadata_and_buffer_bytes':other,'incremental_peak_admission_bytes':usage['allocated_bytes']+needed,'required_group_space_headroom_bytes':needed,'required_inode_headroom':inode_need,'admitted':usage['allocated_bytes']+needed<=cap and quota['block_hard_bytes']-quota['block_used_bytes']-quota['block_in_doubt_bytes']>=needed and quota['files_hard']-quota['files_used']-quota['files_in_doubt']>=inode_need and 22581+usage['file_and_directory_entries']+30*n+2*n+1000<50000}
