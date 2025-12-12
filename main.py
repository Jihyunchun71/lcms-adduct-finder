import sys
import os

# ==========================================
# [안전 장치] 라이브러리 설치 여부 확인
# ==========================================
try:
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    from pyopenms import *
    from molmass import Formula
    from scipy.optimize import curve_fit
except ImportError as e:
    print("\n[!] 필수 라이브러리가 설치되지 않았거나, 경로가 잘못되었습니다.")
    print(f"    에러 상세: {e}")
    print("\n[해결 방법]")
    print("터미널(CMD 또는 PowerShell)에서 아래 명령어를 실행하여 설치해주세요:")
    print(">> pip install pyopenms pandas openpyxl molmass scipy matplotlib numpy")
    print("\n(이미 설치했는데도 이 오류가 뜬다면, 파이썬 실행 경로를 확인해주세요.)")
    sys.exit(1)

# ==========================================
# 1. 환경 설정 및 상수 정의
# ==========================================

# [사용자 설정 필요] mzML 파일들이 들어있는 폴더 경로
# "./"는 현재 이 코드가 있는 폴더를 의미합니다.
# 데이터를 다른 곳에 두셨다면 아래 경로를 수정하세요. (예: r"C:\Data\MyProject")
MZML_DATA_FOLDER = r"./mzml" 

INPUT_EXCEL_FILE = 'file_list.xlsx'  # 입력 엑셀 파일명
INPUT_SHEET_NAME = 'Final'           # 읽어올 시트 이름
EXPORT_PLOT_FOLDER = "EIC_Plots_Export" # 그래프 저장 폴더명
OUTPUT_EXCEL_FILE = "Final_Result_With_Plots.xlsx" # 결과 엑셀 파일명

# 질량 상수
MASS_H = 1.0073     
MASS_Na = 22.9892
MASS_ACN = 41.0265
MASS_FA = 46.0055
MASS_HCOO = 44.9983
MASS_NH4 = 18.0338
MASS_H2O = 18.0106
MASS_E = 0.00054858 # 전자 질량

# Adduct 정의
ADDUCT_DEFINITIONS = {
    # [n=2] Dimer Adducts
    "[2M-2H+Na]-": {"multiplier": 2, "delta": -2 * MASS_H + MASS_Na, "net_charge": -1}, 
    "[2M-H]-":     {"multiplier": 2, "delta": -MASS_H,               "net_charge": -1},
    "[2M+H]+":     {"multiplier": 2, "delta": +MASS_H,               "net_charge": +1},
    # [n=1] Monomer Adducts
    "[M-2H2O+H]+": {"multiplier": 1, "delta": -(2 * MASS_H2O) + MASS_H, "net_charge": +1},
    "[M-3H2O+H]+": {"multiplier": 1, "delta": -(3 * MASS_H2O) + MASS_H, "net_charge": +1},
    "[M-H]-":      {"multiplier": 1, "delta": -MASS_H,                  "net_charge": -1},
    "[M-H2O-H]-":  {"multiplier": 1, "delta": -(MASS_H2O + MASS_H),     "net_charge": -1},
    "[M-H2O+H]+":  {"multiplier": 1, "delta": -MASS_H2O + MASS_H,       "net_charge": +1},
    "[M]-":        {"multiplier": 1, "delta": 0,                        "net_charge": -1}, 
    "[M]+":        {"multiplier": 1, "delta": 0,                        "net_charge": +1},
    "[M+ACN+H]+":  {"multiplier": 1, "delta": +MASS_ACN + MASS_H,       "net_charge": +1},
    "[M+FA-H]-":   {"multiplier": 1, "delta": +(MASS_FA - MASS_H),      "net_charge": -1},
    "[M+H]+":      {"multiplier": 1, "delta": +MASS_H,                  "net_charge": +1},
    "[M+HCOO]-":   {"multiplier": 1, "delta": +MASS_HCOO,               "net_charge": -1},
    "[M+Na]+":     {"multiplier": 1, "delta": +MASS_Na,                 "net_charge": +1},
    "[M+NH4]+":    {"multiplier": 1, "delta": +MASS_NH4,                "net_charge": +1}
}

PPM_TOLERANCE = 10.0

# 그래프 저장 폴더 생성
if not os.path.exists(EXPORT_PLOT_FOLDER):
    os.makedirs(EXPORT_PLOT_FOLDER)

# ==========================================
# 2. 계산 및 분석 함수
# ==========================================

def calculate_target_mz(formula_str, adduct_info):
    try:
        f = Formula(formula_str)
        exact_mass = f.isotope.mass
        mult = adduct_info['multiplier']
        delta = adduct_info['delta']
        charge = adduct_info['net_charge']
        ion_mass = (exact_mass * mult) + delta - (charge * MASS_E)
        return ion_mass / abs(charge)
    except:
        return None

def gaussian_func(x, a, x0, sigma):
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def fit_gaussian_and_score(rt_array, int_array):
    if len(rt_array) < 5 or np.max(int_array) == 0:
        return 0.0, None
    try:
        max_idx = np.argmax(int_array)
        a_guess = int_array[max_idx]
        x0_guess = rt_array[max_idx]
        sigma_guess = 10.0
        
        mask = int_array > (a_guess * 0.1)
        if np.sum(mask) < 4: x_data, y_data = rt_array, int_array
        else: x_data, y_data = rt_array[mask], int_array[mask]

        popt, _ = curve_fit(gaussian_func, x_data, y_data, 
                            p0=[a_guess, x0_guess, sigma_guess], maxfev=2000)
        
        residuals = y_data - gaussian_func(x_data, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y_data - np.mean(y_data))**2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot != 0 else 0.0
        return max(0.0, r_squared), popt
    except:
        return 0.0, None

def save_eic_plot_png(rt_arr, int_arr, fit_params, score, formula, adduct, raw_filename, mz_val):
    try:
        plt.figure(figsize=(6, 4))
        rt_min = rt_arr / 60.0
        plt.plot(rt_min, int_arr, 'b-', label='Raw EIC', linewidth=1.5, alpha=0.7)
        if fit_params is not None:
            fitted_curve = gaussian_func(rt_arr, *fit_params)
            plt.plot(rt_min, fitted_curve, 'r--', label=f'Fit (R2={score:.2f})', linewidth=1.5)

        plt.title(f"{formula} {adduct}\nFile: {raw_filename}", fontsize=10)
        plt.xlabel("Retention Time (min)")
        plt.ylabel("Intensity")
        plt.legend(loc='upper right', fontsize='small')
        plt.grid(True, linestyle=':', alpha=0.6)
        plt.tight_layout()
        
        safe_form = "".join(c for c in formula if c.isalnum())
        safe_add = adduct.replace("[", "").replace("]", "").replace("+", "p").replace("-", "m")
        save_name = f"{safe_form}_{safe_add}_{raw_filename}.png"
        
        save_path = os.path.join(EXPORT_PLOT_FOLDER, save_name)
        plt.savefig(save_path, dpi=100)
    except Exception as e:
        print(f"Plot error: {e}")
    finally:
        plt.close()

def process_ms_file(file_path, target_list, output_data, mode):
    filename = os.path.basename(file_path)
    print(f"Processing: {filename} ({mode})")
    
    exp = MSExperiment()
    try:
        MzMLFile().load(file_path, exp)
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        return

    ms1_data = [] 
    ms2_precursors = []
    for spec in exp:
        if spec.getMSLevel() == 1:
            ms1_data.append((spec.getRT(), spec.get_peaks()[0], spec.get_peaks()[1]))
        elif spec.getMSLevel() == 2:
            precursors = spec.getPrecursors()
            if precursors: ms2_precursors.append(precursors[0].getMZ())

    for idx, row in target_list.iterrows():
        formula = row['Formula']
        for adduct_name, adduct_info in ADDUCT_DEFINITIONS.items():
            if (mode == 'POS' and adduct_info['net_charge'] < 0) or \
               (mode == 'NEG' and adduct_info['net_charge'] > 0):
                continue
                
            target_mz = calculate_target_mz(formula, adduct_info)
            if target_mz is None: continue

            eic_rt = []
            eic_int = []
            ppm_error = target_mz * (PPM_TOLERANCE / 1e6)
            
            for rt, mz_arr, int_arr in ms1_data:
                mask = (mz_arr >= target_mz - ppm_error) & (mz_arr <= target_mz + ppm_error)
                if np.any(mask):
                    eic_rt.append(rt)
                    eic_int.append(np.sum(int_arr[mask]))
                else:
                    eic_rt.append(rt)
                    eic_int.append(0.0)
            
            eic_rt = np.array(eic_rt)
            eic_int = np.array(eic_int)
            max_intensity = np.max(eic_int)
            total_area = np.sum(eic_int)
            
            best_rt_min = 0.0
            gauss_score = 0.0
            quality_label = "Noise"
            
            if max_intensity > 1000:
                best_rt_min = eic_rt[np.argmax(eic_int)] / 60.0
                gauss_score, fit_params = fit_gaussian_and_score(eic_rt, eic_int)
                save_eic_plot_png(eic_rt, eic_int, fit_params, gauss_score, 
                                  formula, adduct_name, filename, target_mz)
                
                if gauss_score > 0.8: quality_label = "Excellent"
                elif gauss_score > 0.5: quality_label = "Good"
                else: quality_label = "Poor Shape"

            has_ms2 = False
            ppm_error_ms2 = target_mz * (PPM_TOLERANCE / 1e6)
            for prec_mz in ms2_precursors:
                if abs(prec_mz - target_mz) <= ppm_error_ms2:
                    has_ms2 = True
                    break

            output_data.append({
                'RawFile': filename,
                'Mode': mode,
                'Formula': formula,
                'Adduct': adduct_name,
                'mz_theoretical': target_mz,
                'RT_min': round(best_rt_min, 3),
                'Intensity': max_intensity,
                'Area': total_area,
                'GaussianScore': round(gauss_score, 3),
                'PeakQuality': quality_label,
                'HasMS2': has_ms2
            })

# ==========================================
# 3. 메인 실행 (폴더 경로 적용)
# ==========================================

def main():
    if not os.path.exists(INPUT_EXCEL_FILE):
        print(f"오류: 엑셀 파일 '{INPUT_EXCEL_FILE}'이 스크립트와 같은 폴더에 없습니다.")
        return

    print("엑셀 파일을 읽는 중...")
    try:
        meta_data = pd.read_excel(INPUT_EXCEL_FILE, sheet_name=INPUT_SHEET_NAME)
    except Exception as e:
        print(f"오류: {e}")
        return

    all_results = []

    for idx, row in meta_data.iterrows():
        raw_filename = str(row['RawFile']).strip() # 엑셀에 적힌 파일명
        mode = row['Mode']
        
        # [자동 보정] 확장자가 .raw로 되어있으면 .mzML로 바꿔서 찾기 시도
        if raw_filename.lower().endswith('.raw'):
            mzml_filename = raw_filename[:-4] + ".mzML"
        elif not raw_filename.lower().endswith('.mzML'):
            mzml_filename = raw_filename + ".mzML"
        else:
            mzml_filename = raw_filename
            
        # [경로 결합] 지정된 폴더 + 파일명
        full_file_path = os.path.join(MZML_DATA_FOLDER, mzml_filename)
        
        # 화학식 파싱
        raw_formulas = row.iloc[2:].values
        formulas = [str(f).strip() for f in raw_formulas if pd.notna(f) and str(f).strip() != '']
        
        print(f"\n[{idx+1}/{len(meta_data)}] File: {mzml_filename}")
        
        if os.path.exists(full_file_path):
            target_df = pd.DataFrame({'Formula': formulas})
            process_ms_file(full_file_path, target_df, all_results, mode)
        else:
            print(f"  - (Error) 파일을 찾을 수 없습니다.")
            print(f"  - 경로 확인: {full_file_path}")

    if all_results:
        print("\n데이터 저장 중...")
        df_results = pd.DataFrame(all_results)
        
        with pd.ExcelWriter(OUTPUT_EXCEL_FILE, engine='openpyxl') as writer:
            # 1. 전체 Raw Data 시트
            df_results.to_excel(writer, sheet_name='All_Features', index=False)
            
            unique_formulas = df_results['Formula'].unique()
            for formula in unique_formulas:
                f_data = df_results[df_results['Formula'] == formula]
                
                # 피벗 테이블 생성
                pivot_area = f_data.pivot_table(index='RawFile', columns='Adduct', values='Area')
                pivot_rt = f_data.pivot_table(index='RawFile', columns='Adduct', values='RT_min')
                
                # 시트명 생성
                safe_name = "".join(c for c in formula if c.isalnum())[:30]
                
                # -------------------------------------------------------
                # [1] Area Table 저장
                # -------------------------------------------------------
                pivot_area.to_excel(writer, sheet_name=safe_name, startrow=0)
                writer.sheets[safe_name].cell(row=1, column=1).value = "Area Table"
                
                # -------------------------------------------------------
                # [2] Retention Time Table 저장
                # -------------------------------------------------------
                
                # 다음 표가 시작될 위치 계산: (Area 표 길이) + (여백 3줄)
                current_row = len(pivot_area) + 3
                
                writer.sheets[safe_name].cell(row=current_row, column=1).value = "Retention Time (min)"
                pivot_rt.to_excel(writer, sheet_name=safe_name, startrow=current_row+1)

        print(f"\n[완료] 결과 엑셀: {OUTPUT_EXCEL_FILE}")
    else:
        print("\n[알림] 저장할 데이터가 없습니다.")

if __name__ == "__main__":
    main()