import pandas as pd
import pytest

from biogroups import compute_biological_groups


def write_metadata(tmp_path):
    extractmetadata = tmp_path / 'extractmetadata.csv'
    pd.DataFrame({
        'Sample_Code': ['Blank', 'MB1', 'MB2', 'MB3'],
        'Biological_Group': ['Blanks', '0um_Ce', '0um_Ce', '250um_Ce'],
    }).to_csv(extractmetadata, index=False)

    samplelist = tmp_path / 'samplelist.csv'
    pd.DataFrame({
        'Injection': ['inj_blank1', 'inj_mb1_r1', 'inj_mb2_r1', 'inj_mb3_r1'],
        'Sample_Code': ['Blank', 'MB1', 'MB2', 'MB3'],
    }).to_csv(samplelist, index=False)

    return extractmetadata, samplelist


def write_peaktable(tmp_path, injections):
    # Canonical 3-header-row peak table: row2 holds injection names.
    path = tmp_path / 'peaktable.csv'
    with open(path, 'w') as f:
        f.write(',,,' + ','.join('' for _ in injections) + '\n')
        f.write(',,,' + ','.join('' for _ in injections) + '\n')
        f.write('Compound,m/z,Retention time (min),' + ','.join(injections) + '\n')
        f.write('c1,100.0,1.0,' + ','.join('1000' for _ in injections) + '\n')
    return path


def test_returns_distinct_groups_in_first_seen_order(tmp_path):
    extractmetadata, samplelist = write_metadata(tmp_path)
    peaktable = write_peaktable(tmp_path, ['inj_blank1', 'inj_mb1_r1', 'inj_mb2_r1', 'inj_mb3_r1'])

    groups, unresolved = compute_biological_groups(extractmetadata, samplelist, peaktable)

    assert groups == ['Blanks', '0um_Ce', '250um_Ce']
    assert unresolved == []


def test_injection_missing_from_metadata_is_reported_not_raised(tmp_path):
    extractmetadata, samplelist = write_metadata(tmp_path)
    peaktable = write_peaktable(tmp_path, ['inj_blank1', 'inj_unknown_sample'])

    groups, unresolved = compute_biological_groups(extractmetadata, samplelist, peaktable)

    assert groups == ['Blanks']
    assert unresolved == ['inj_unknown_sample']


def test_duplicate_injection_row_takes_first_match(tmp_path):
    extractmetadata, samplelist = write_metadata(tmp_path)
    # Append a second, duplicate Sample_Code/Injection row -- makes the
    # joined metadata's .loc[] return a Series instead of a scalar for that
    # injection.
    df = pd.read_csv(samplelist)
    df = pd.concat([df, pd.DataFrame({'Injection': ['inj_blank1'], 'Sample_Code': ['Blank']})], ignore_index=True)
    df.to_csv(samplelist, index=False)

    peaktable = write_peaktable(tmp_path, ['inj_blank1'])

    groups, unresolved = compute_biological_groups(extractmetadata, samplelist, peaktable)

    assert groups == ['Blanks']
    assert unresolved == []


def test_bad_metadata_join_raises_value_error(tmp_path):
    extractmetadata = tmp_path / 'extractmetadata.csv'
    pd.DataFrame({'NotSampleCode': ['x']}).to_csv(extractmetadata, index=False)
    samplelist = tmp_path / 'samplelist.csv'
    pd.DataFrame({'Injection': ['x'], 'Sample_Code': ['x']}).to_csv(samplelist, index=False)
    peaktable = write_peaktable(tmp_path, ['x'])

    with pytest.raises(ValueError):
        compute_biological_groups(extractmetadata, samplelist, peaktable)
