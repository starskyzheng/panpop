fn main() {
    let args = std::env::args();
    println!("{:?}", args);
    let inrgfa_file:String = String::from_str(&args[1]).expect("Error input rgfa file");
    let outchrs_file:String = String::from_str(&args[2]).expect("Error outchrs_file file");
    println!("{}", inrgfa_file);
}