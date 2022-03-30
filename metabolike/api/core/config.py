from pydantic import BaseSettings, SecretStr


class DatabaseConfig(BaseSettings):
    metabolike_db_uri: str
    metabolike_db_user: str
    metabolike_db_password: SecretStr
    metabolike_db_name: str

    class Config:
        env_file = "metabolike.env"


db_conf = DatabaseConfig()
